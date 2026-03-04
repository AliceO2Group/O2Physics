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

/// \file CorrSparse.cxx
/// \brief Provides a sparse with usefull two particle correlation info
/// \author Thor Jensen (thor.kjaersgaard.jensen@cern.ch)

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWCumulant.h"
#include "PWGCF/GenericFramework/Core/GFWPowerArray.h"
#include "PWGCF/GenericFramework/Core/GFWWeights.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "FT0Base/Geometry.h"
#include "FV0Base/Geometry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>

#include "TRandom3.h"

#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
namespace o2::aod
{
namespace corrsparse
{
DECLARE_SOA_COLUMN(Multiplicity, multiplicity, int);
}
DECLARE_SOA_TABLE(Multiplicity, "AOD", "MULTIPLICITY",
                  corrsparse::Multiplicity);

} // namespace o2::aod

// define the filtered collisions and tracks
#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

// static constexpr float LongArrayFloat[3][20] = {{1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {2.1, 2.2, 2.3, -2.1, -2.2, -2.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {3.1, 3.2, 3.3, -3.1, -3.2, -3.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}};
// static constexpr int LongArrayInt[3][20] = {{1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1}, {2, 2, 2, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1}, {3, 3, 3, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1}};

struct CorrSparse {
  Service<ccdb::BasicCCDBManager> ccdb;

  O2_DEFINE_CONFIGURABLE(cfgUseAdditionalEventCut, bool, false, "Use additional event cut on mult correlations")
  O2_DEFINE_CONFIGURABLE(cfgZVtxCut, float, 10.0f, "Accepted z-vertex range")

  struct : ConfigurableGroup{
             O2_DEFINE_CONFIGURABLE(cfgPtCutMin, float, 0.2f, "minimum accepted track pT")
               O2_DEFINE_CONFIGURABLE(cfgPtCutMax, float, 10.0f, "maximum accepted track pT")
                 O2_DEFINE_CONFIGURABLE(cfgEtaCut, float, 0.8f, "Eta cut")
                   O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5f, "max chi2 per TPC clusters")
                     O2_DEFINE_CONFIGURABLE(cfgCutTPCclu, float, 50.0f, "minimum TPC clusters")
                       O2_DEFINE_CONFIGURABLE(cfgCutTPCCrossedRows, float, 70.0f, "minimum TPC crossed rows")
                         O2_DEFINE_CONFIGURABLE(cfgCutITSclu, float, 5.0f, "minimum ITS clusters")
                           O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "max DCA to vertex z")

           } cfgTrackCuts;

  struct : ConfigurableGroup{

             O2_DEFINE_CONFIGURABLE(processFV0, bool, true, "Process FV0 correlations")
               O2_DEFINE_CONFIGURABLE(processFT0A, bool, true, "Process FT0A correlations")
                 O2_DEFINE_CONFIGURABLE(processFT0C, bool, true, "Process FT0C correlations")
                   O2_DEFINE_CONFIGURABLE(processMFT, bool, true, "Process MFT correlations")

           } cfgDetectorConfig;

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgPtCutMinMFT, float, 0.5f, "minimum accepted MFT track pT")
    O2_DEFINE_CONFIGURABLE(cfgPtCutMaxMFT, float, 10.0f, "maximum accepted MFT track pT")
    O2_DEFINE_CONFIGURABLE(etaMftTrackMin, float, -3.6, "Minimum eta for MFT track")
    O2_DEFINE_CONFIGURABLE(etaMftTrackMax, float, -2.5, "Maximum eta for MFT track")
    O2_DEFINE_CONFIGURABLE(nClustersMftTrack, int, 5, "Minimum number of clusters for MFT track")
    Configurable<int> cutBestCollisionId{"cutBestCollisionId", 0, "cut on the best collision Id used in a filter"};
    Configurable<float> etaMftTrackMaxFilter{"etaMftTrackMaxFilter", -2.0f, "Maximum value for the eta of MFT tracks when used in filter"};
    Configurable<float> etaMftTrackMinFilter{"etaMftTrackMinFilter", -3.9f, "Minimum value for the eta of MFT tracks when used in filter"};
    Configurable<float> mftMaxDCAxy{"mftMaxDCAxy", 2.0f, "Cut on dcaXY for MFT tracks"};
    Configurable<float> mftMaxDCAz{"mftMaxDCAz", 2.0f, "Cut on dcaZ for MFT tracks"};
  } cfgMftConfig;

  struct : ConfigurableGroup{
             O2_DEFINE_CONFIGURABLE(cfgMinMult, int, 0, "Minimum multiplicity for collision")
               O2_DEFINE_CONFIGURABLE(cfgMaxMult, int, 10, "Maximum multiplicity for collision")
                 O2_DEFINE_CONFIGURABLE(cfgEvSelkNoSameBunchPileup, bool, false, "rejects collisions which are associated with the same found-by-T0 bunch crossing")
                   O2_DEFINE_CONFIGURABLE(cfgEvSelkNoITSROFrameBorder, bool, false, "reject events at ITS ROF border")
                     O2_DEFINE_CONFIGURABLE(cfgEvSelkNoTimeFrameBorder, bool, false, "reject events at TF border")
                       O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodZvtxFT0vsPV, bool, false, "removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference, use this cut at low multiplicities with caution")
                         O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInTimeRangeStandard, bool, false, "no collisions in specified time range")
                           O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodITSLayer0123, bool, true, "cut time intervals with dead ITS layers 0,1,2,3")
                             O2_DEFINE_CONFIGURABLE(cfgEvSelkIsGoodITSLayersAll, bool, true, "cut time intervals with dead ITS staves")
                               O2_DEFINE_CONFIGURABLE(cfgEvSelkNoCollInRofStandard, bool, false, "no other collisions in this Readout Frame with per-collision multiplicity above threshold")
                                 O2_DEFINE_CONFIGURABLE(cfgEvSelkNoHighMultCollInPrevRof, bool, false, "veto an event if FT0C amplitude in previous ITS ROF is above threshold")
                                   O2_DEFINE_CONFIGURABLE(cfgEvSelMultCorrelation, bool, true, "Multiplicity correlation cut")
                                     O2_DEFINE_CONFIGURABLE(cfgEvSelV0AT0ACut, bool, true, "V0A T0A 5 sigma cut")
                                       O2_DEFINE_CONFIGURABLE(cfgEvSelOccupancy, bool, true, "Occupancy cut")
                                         O2_DEFINE_CONFIGURABLE(cfgCutOccupancyHigh, int, 2000, "High cut on TPC occupancy")
                                           O2_DEFINE_CONFIGURABLE(cfgCutOccupancyLow, int, 0, "Low cut on TPC occupancy")

           } cfgEventSelection;

  struct : ConfigurableGroup{

             O2_DEFINE_CONFIGURABLE(cfgRejectFT0AInside, bool, false, "Rejection of inner ring channels of the FT0A detector")
               O2_DEFINE_CONFIGURABLE(cfgRejectFT0AOutside, bool, false, "Rejection of outer ring channels of the FT0A detector")
                 O2_DEFINE_CONFIGURABLE(cfgRejectFT0CInside, bool, false, "Rejection of inner ring channels of the FT0C detector")
                   O2_DEFINE_CONFIGURABLE(cfgRejectFT0COutside, bool, false, "Rejection of outer ring channels of the FT0C detector")
                     O2_DEFINE_CONFIGURABLE(cfgRemapFT0ADeadChannels, bool, false, "If true, remap FT0A channels 60-63 to amplitudes from 92-95 respectively")
                       O2_DEFINE_CONFIGURABLE(cfgRemapFT0CDeadChannels, bool, false, "If true, remap FT0C channels 177->145, 176->144, 178->146, 179->147, 139->115")

           } cfgFITConfig;

  O2_DEFINE_CONFIGURABLE(cfgMinMixEventNum, int, 5, "Minimum number of events to mix")
  O2_DEFINE_CONFIGURABLE(cfgMergingCut, float, 0.02, "Merging cut on track merge")
  O2_DEFINE_CONFIGURABLE(cfgApplyTwoTrackEfficiency, bool, true, "Apply two track efficiency for tpc tpc")
  O2_DEFINE_CONFIGURABLE(cfgRadiusLow, float, 0.8, "Low radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgRadiusHigh, float, 2.5, "High radius for merging cut")
  O2_DEFINE_CONFIGURABLE(cfgSampleSize, double, 10, "Sample size for mixed event")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgCentralityWeight, std::string, "", "CCDB path to centrality weight object")
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, bool, false, "Use local efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgUseEventWeights, bool, false, "Use event weights for mixed event")

  struct : ConfigurableGroup {
    O2_DEFINE_CONFIGURABLE(cfgMultCentHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 10.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultCentLowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultT0CCutEnabled, bool, false, "Enable Global multiplicity vs T0C centrality cut")
    Configurable<std::vector<double>> cfgMultT0CCutPars{"cfgMultT0CCutPars", std::vector<double>{143.04, -4.58368, 0.0766055, -0.000727796, 2.86153e-06, 23.3108, -0.36304, 0.00437706, -4.717e-05, 1.98332e-07}, "Global multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultPVT0CCutEnabled, bool, false, "Enable PV multiplicity vs T0C centrality cut")
    Configurable<std::vector<double>> cfgMultPVT0CCutPars{"cfgMultPVT0CCutPars", std::vector<double>{195.357, -6.15194, 0.101313, -0.000955828, 3.74793e-06, 30.0326, -0.43322, 0.00476265, -5.11206e-05, 2.13613e-07}, "PV multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultMultPVHighCutFunction, std::string, "[0]+[1]*x + 5.*([2]+[3]*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultPVLowCutFunction, std::string, "[0]+[1]*x - 5.*([2]+[3]*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultGlobalPVCutEnabled, bool, false, "Enable global multiplicity vs PV multiplicity cut")
    Configurable<std::vector<double>> cfgMultGlobalPVCutPars{"cfgMultGlobalPVCutPars", std::vector<double>{-0.140809, 0.734344, 2.77495, 0.0165935}, "PV multiplicity vs T0C centrality cut parameter values"};
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0AHighCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + 4.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0ALowCutFunction, std::string, "[0] + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x - 3.*([5] + [6]*x + [7]*x*x + [8]*x*x*x + [9]*x*x*x*x)", "Functional for multiplicity correlation cut");
    O2_DEFINE_CONFIGURABLE(cfgMultMultV0ACutEnabled, bool, false, "Enable global multiplicity vs V0A multiplicity cut")
    Configurable<std::vector<double>> cfgMultMultV0ACutPars{"cfgMultMultV0ACutPars", std::vector<double>{534.893, 184.344, 0.423539, -0.00331436, 5.34622e-06, 871.239, 53.3735, -0.203528, 0.000122758, 5.41027e-07}, "Global multiplicity vs V0A multiplicity cut parameter values"};
    std::vector<double> multT0CCutPars;
    std::vector<double> multPVT0CCutPars;
    std::vector<double> multGlobalPVCutPars;
    std::vector<double> multMultV0ACutPars;
    TF1* fMultPVT0CCutLow = nullptr;
    TF1* fMultPVT0CCutHigh = nullptr;
    TF1* fMultT0CCutLow = nullptr;
    TF1* fMultT0CCutHigh = nullptr;
    TF1* fMultGlobalPVCutLow = nullptr;
    TF1* fMultGlobalPVCutHigh = nullptr;
    TF1* fMultMultV0ACutLow = nullptr;
    TF1* fMultMultV0ACutHigh = nullptr;
    TF1* fT0AV0AMean = nullptr;
    TF1* fT0AV0ASigma = nullptr;
  } cfgFuncParas;

  SliceCache cache;
  SliceCache cacheNch;

  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {72, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt axis for histograms"};
  ConfigurableAxis axisAmbiguity{"axisAmbiguity", {100, 0, 100}, "MFT track ambiguity axis for histograms"};

  ConfigurableAxis axisEtaMft{"axisEtaMFT", {40, -3.6, -2.4}, "eta axis for MFT tracks in histograms"};
  ConfigurableAxis axisEtaFt0a{"axisEtaFT0A", {40, 3.5, 4.9}, "eta axis for FT0A in histograms"};
  ConfigurableAxis axisEtaFt0c{"axisEtaFT0C", {40, -3.3, -2.1}, "eta axis for FT0C in histograms"};
  ConfigurableAxis axisEtaFv0{"axisEtaFV0", {40, 2.2, 5.1}, "eta axis for FV0 in histograms"};

  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {48, -1.6, 1.6}, "delta eta axis for histograms"};
  ConfigurableAxis axisDeltaEtaTpcMft{"axisDeltaEtaTPCMFT", {48, 1.6, 4.6}, "delta eta axis for TPC-MFT histograms"};
  ConfigurableAxis axisDeltaEtaTpcFv0{"axisDeltaEtaTPCFV0", {48, -1.7, -5.9}, "delta eta axis for TPC-FV0 histograms"};
  ConfigurableAxis axisDeltaEtaTpcFt0a{"axisDeltaEtaTPCFT0A", {48, -5.7, -2.7}, "delta eta axis for TPC-FT0A histograms"};
  ConfigurableAxis axisDeltaEtaTpcFt0c{"axisDeltaEtaTPCFT0C", {48, 1.3, 4.1}, "delta eta axis for TPC-FT0C histograms"};
  ConfigurableAxis axisDeltaEtaMftFt0c{"axisDeltaEtaMFTFT0C", {48, -2.0, 0.6}, "delta eta axis for MFT-FT0C histograms"};
  ConfigurableAxis axisDeltaEtaMftFt0a{"axisDeltaEtaMFTFT0A", {48, -8.5, -5.9}, "delta eta axis for MFT-FT0A histograms"};
  ConfigurableAxis axisDeltaEtaMftFv0{"axisDeltaEtaMFTFV0", {48, -8.6, -4.7}, "delta eta axis for MFT-FV0 histograms"};

  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity / centrality axis for histograms"};
  ConfigurableAxis vtxMix{"vtxMix", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "vertex axis for mixed event histograms"};
  ConfigurableAxis multMix{"multMix", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity / centrality axis for mixed event histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};

  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisNsigmaITS{"axisNsigmaITS", {80, -5, 5}, "nsigmaITS axis"};
  ConfigurableAxis axisTpcSignal{"axisTpcSignal", {250, 0, 250}, "dEdx axis for TPC"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  ConfigurableAxis axisAmplitudeFt0a{"axisAmplitudeFt0a", {5000, 0, 1000}, "FT0A amplitude"};
  ConfigurableAxis axisChannelFt0aAxis{"axisChannelFt0aAxis", {96, 0.0, 96.0}, "FT0A channel"};

  Configurable<std::string> cfgGainEqPath{"cfgGainEqPath", "Analysis/EventPlane/GainEq", "CCDB path for gain equalization constants"};
  Configurable<int> cfgCorrLevel{"cfgCorrLevel", 1, "calibration step: 0 = no corr, 1 = gain corr"};
  ConfigurableAxis cfgaxisFITamp{"cfgaxisFITamp", {1000, 0, 5000}, ""};
  AxisSpec axisFit{cfgaxisFITamp, "fit amplitude"};
  AxisSpec axisChID = {220, 0, 220};

  // make the filters and cuts.
  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgZVtxCut);
  Filter trackFilter = (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaCut) && (cfgTrackCuts.cfgPtCutMin < aod::track::pt) && (cfgTrackCuts.cfgPtCutMax > aod::track::pt) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgTrackCuts.cfgCutChi2prTPCcls) && (aod::track::dcaZ < cfgTrackCuts.cfgCutDCAz);

  Filter mftTrackEtaFilter = ((aod::fwdtrack::eta < cfgMftConfig.etaMftTrackMaxFilter) && (aod::fwdtrack::eta > cfgMftConfig.etaMftTrackMinFilter));

  // Filters below will be used for uncertainties
  Filter mftTrackCollisionIdFilter = (aod::fwdtrack::bestCollisionId >= 0);
  Filter mftTrackDcaXYFilter = (nabs(aod::fwdtrack::bestDCAXY) < cfgMftConfig.mftMaxDCAxy);
  // Filter mftTrackDcaZFilter = (nabs(aod::fwdtrack::bestDCAZ) < cfgMftConfig.mftMaxDCAz);

  // using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel>>; // aod::CentFT0Cs
  // using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra>>;

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSel, aod::Mults>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA>>;
  using FilteredMftTracks = soa::Filtered<aod::MFTTracks>;
  using Reassociated2DMftTracks = aod::BestCollisionsFwd;

  Preslice<AodTracks> perColGlobal = aod::track::collisionId;

  // FT0 geometry
  o2::ft0::Geometry ft0Det;
  o2::fv0::Geometry* fv0Det{};
  static constexpr uint64_t Ft0IndexA = 96;
  std::vector<o2::detectors::AlignParam>* offsetFT0;
  std::vector<o2::detectors::AlignParam>* offsetFV0;

  std::vector<float> cstFT0RelGain{};

  // Corrections
  TH3D* mEfficiency = nullptr;
  TH1D* mCentralityWeight = nullptr;
  bool correctionsLoaded = false;

  // Define the outputs

  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};

  HistogramRegistry registry{"registry"};

  // Define Global variables

  enum EventCutTypes {
    kFilteredEvents = 0,
    kAfterSel8,
    kUseNoTimeFrameBorder,
    kUseNoITSROFrameBorder,
    kUseNoSameBunchPileup,
    kUseGoodZvtxFT0vsPV,
    kUseNoCollInTimeRangeStandard,
    kUseGoodITSLayersAll,
    kUseGoodITSLayer0123,
    kUseNoCollInRofStandard,
    kUseNoHighMultCollInPrevRof,
    kUseOccupancy,
    kUseMultCorrCut,
    kUseT0AV0ACut,
    kUseVertexITSTPC,
    kUseTVXinTRD,
    kNEventCuts
  };

  enum MftTrackAmbiguityStep {
    AllMftTracks = 0,
    AfterTrackSelection,
    NumberOfAmbiguousTracks,
    NumberOfNonAmbiguousTracks,
    NMftAmbiguitySteps
  };

  enum ReassociationMftTracks {
    NotReassociatedMftTracks = 0,
    ReassociatedMftTracks,
    NReassociationMftTracksSteps
  };

  enum EventType {
    SameEvent = 1,
    MixedEvent = 3
  };

  enum FITIndex {
    kFT0A = 0,
    kFT0C,
    kFV0
  };

  enum DetectorChannels {
    kFT0AInnerRingMin = 0,
    kFT0AInnerRingMax = 31,
    kFT0AOuterRingMin = 32,
    kFT0AOuterRingMax = 95,
    kFT0CInnerRingMin = 96,
    kFT0CInnerRingMax = 143,
    kFT0COuterRingMin = 144,
    kFT0COuterRingMax = 207
  };

  std::array<std::array<int, 1>, 16> eventCuts;

  void init(InitContext&)
  {

    const AxisSpec axisPhi{72, 0.0, constants::math::TwoPI, "#varphi"};
    const AxisSpec axisEta{40, -1., 1., "#eta"};
    const AxisSpec axisEtaFull{90, -4., 5., "#eta"};

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
    fv0Det = o2::fv0::Geometry::instance(o2::fv0::Geometry::eUninitialized);

    LOGF(info, "Starting init");

    // Event Counter
    if ((doprocessSameTpcFIT || doprocessSameTpcMft || doprocessSameTPC || doprocessSameMFTFIT || doprocessSameTpcMftReassociated2D || doprocessSameTpcMftReassociated3D) && cfgUseAdditionalEventCut) {
      registry.add("hEventCountSpecific", "Number of Event;; Count", {HistType::kTH1D, {{13, 0, 13}}});
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(1, "after sel8");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(2, "kNoSameBunchPileup");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(3, "kNoITSROFrameBorder");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(4, "kNoTimeFrameBorder");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(5, "kIsGoodZvtxFT0vsPV");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(6, "kNoCollInTimeRangeStandard");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(7, "kIsGoodITSLayer0123");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(8, "kIsGoodITSLayersAll");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(9, "kNoCollInRofStandard");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(10, "kNoHighMultCollInPrevRof");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(11, "occupancy");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(12, "MultCorrelation");
      registry.get<TH1>(HIST("hEventCountSpecific"))->GetXaxis()->SetBinLabel(13, "cfgEvSelV0AT0ACut");
    }

    if (doprocessSameTpcMftReassociated2D || doprocessSameTpcMftReassociated3D) {
      registry.add("hEventCountMftReassoc", "Number of Event;; Count", {HistType::kTH1D, {{5, 0, 5}}});
      registry.get<TH1>(HIST("hEventCountMftReassoc"))->GetXaxis()->SetBinLabel(1, "all MFT tracks");
      registry.get<TH1>(HIST("hEventCountMftReassoc"))->GetXaxis()->SetBinLabel(2, "MFT tracks after selection");
      registry.get<TH1>(HIST("hEventCountMftReassoc"))->GetXaxis()->SetBinLabel(3, "ambiguous MFT tracks");
      registry.get<TH1>(HIST("hEventCountMftReassoc"))->GetXaxis()->SetBinLabel(4, "non-ambiguous MFT tracks");
      registry.get<TH1>(HIST("hEventCountMftReassoc"))->GetXaxis()->SetBinLabel(5, "Reassociated MFT tracks");

      registry.add("ReassociatedMftTracks", "Reassociated MFT tracks", {HistType::kTH1D, {{3, 0, 3}}});
      registry.get<TH1>(HIST("ReassociatedMftTracks"))->GetXaxis()->SetBinLabel(1, "Not Reassociated MFT tracks");
      registry.get<TH1>(HIST("ReassociatedMftTracks"))->GetXaxis()->SetBinLabel(2, "Reassociated MFT tracks");
    }

    // Make histograms to check the distributions after cuts
    if (doprocessSameTpcFIT || doprocessSameTpcMft || doprocessSameTPC || doprocessSameMFTFIT || doprocessSameTpcMftReassociated2D || doprocessSameTpcMftReassociated3D) {
      registry.add("Phi", "Phi", {HistType::kTH1D, {axisPhi}});
      if (doprocessSameMFTFIT) {
        registry.add("Eta", "EtaMFT", {HistType::kTH1D, {axisEtaMft}});
      }
      if (doprocessSameTpcFIT || doprocessSameTPC || doprocessSameTpcMftReassociated2D || doprocessSameTpcMftReassociated3D) {
        registry.add("Eta", "Eta", {HistType::kTH1D, {axisEta}});
      }
      registry.add("EtaCorrected", "EtaCorrected", {HistType::kTH1D, {axisEta}});
      registry.add("pT", "pT", {HistType::kTH1D, {axisPtTrigger}});
      registry.add("pTCorrected", "pTCorrected", {HistType::kTH1D, {axisPtTrigger}});
      registry.add("Nch", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}});
      registry.add("Nch_used", "N_{ch}", {HistType::kTH1D, {axisMultiplicity}}); // histogram to see how many events are in the same and mixed event
      registry.add("zVtx", "zVtx", {HistType::kTH1D, {axisVertex}});
      registry.add("zVtx_used", "zVtx_used", {HistType::kTH1D, {axisVertex}});

      if (doprocessSameTpcFIT || doprocessSameMFTFIT) {
        registry.add("FT0Amp", "", {HistType::kTH2F, {axisChID, axisFit}});
        registry.add("FV0Amp", "", {HistType::kTH2F, {axisChID, axisFit}});
        registry.add("FT0AmpCorrect", "", {HistType::kTH2F, {axisChID, axisFit}});
        registry.add("EtaPhi", "", {HistType::kTH2F, {axisEtaFull, axisPhi}});
      }
    }

    if (doprocessSameTpcFIT) {

      if (cfgDetectorConfig.processFT0A) {
        registry.add("deltaEta_deltaPhi_same_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}}); // check to see the delta eta and delta phi distribution
        registry.add("deltaEta_deltaPhi_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0a}});
        registry.add("Assoc_amp_same_TPC_FT0A", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
        registry.add("Assoc_amp_mixed_TPC_FT0A", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
        registry.add("Trig_hist_TPC_FT0A", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
      }
      if (cfgDetectorConfig.processFT0C) {
        registry.add("deltaEta_deltaPhi_same_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}}); // check to see the delta eta and delta phi distribution
        registry.add("deltaEta_deltaPhi_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFt0c}});
        registry.add("Assoc_amp_same_TPC_FT0C", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
        registry.add("Assoc_amp_mixed_TPC_FT0C", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
        registry.add("Trig_hist_TPC_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
      }
      if (cfgDetectorConfig.processFV0) {
        registry.add("deltaEta_deltaPhi_same_TPC_FV0", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFv0}}); // check to see the delta eta and delta phi distribution
        registry.add("deltaEta_deltaPhi_same_TPC_FV0", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFv0}}); // check to see the delta eta and delta phi distribution
        registry.add("deltaEta_deltaPhi_mixed_TPC_FV0", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcFv0}});
        registry.add("Trig_hist_FT0A_FT0C", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
      }
    }

    if (doprocessSameMFTFIT) {

      if (cfgDetectorConfig.processFT0A) {
        registry.add("deltaEta_deltaPhi_same_MFT_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaMftFt0a}}); // check to see the delta eta and delta phi distribution
        registry.add("deltaEta_deltaPhi_mixed_MFT_FT0A", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaMftFt0a}});
        registry.add("Assoc_amp_same", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
        registry.add("Assoc_amp_mixed", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
        registry.add("Trig_hist", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
      }
      if (cfgDetectorConfig.processFT0C) {
        registry.add("deltaEta_deltaPhi_same_MFT_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaMftFt0c}}); // check to see the delta eta and delta phi distribution
        registry.add("deltaEta_deltaPhi_mixed_MFT_FT0C", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaMftFt0c}});
        registry.add("Assoc_amp_same", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
        registry.add("Assoc_amp_mixed", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
        registry.add("Trig_hist", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
      }
      if (cfgDetectorConfig.processFV0) {
        registry.add("deltaEta_deltaPhi_same_MFT_FV0", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaMftFv0}}); // check to see the delta eta and delta phi distribution
        registry.add("deltaEta_deltaPhi_mixed_MFT_FV0", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaMftFv0}});
        registry.add("Trig_hist", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
        registry.add("Assoc_amp_same", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
        registry.add("Assoc_amp_mixed", "", {HistType::kTH2D, {axisChannelFt0aAxis, axisAmplitudeFt0a}});
      }
    }

    if (doprocessSameTpcMft || doprocessSameTpcMftReassociated2D || doprocessSameTpcMftReassociated3D) {
      registry.add("deltaEta_deltaPhi_same_TPC_MFT", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcMft}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_TPC_MFT", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEtaTpcMft}});
      registry.add("Trig_hist_TPC_MFT", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
    }

    if (doprocessSameTPC) {
      registry.add("deltaEta_deltaPhi_same_TPC_TPC", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}}); // check to see the delta eta and delta phi distribution
      registry.add("deltaEta_deltaPhi_mixed_TPC_TPC", "", {HistType::kTH2D, {axisDeltaPhi, axisDeltaEta}});
      registry.add("Trig_hist_TPC_TPC", "", {HistType::kTHnSparseF, {{axisSample, axisVertex, axisPtTrigger}}});
    }

    registry.add("eventcount", "bin", {HistType::kTH1F, {{4, 0, 4, "bin"}}}); // histogram to see how many events are in the same and mixed event

    LOGF(info, "Initializing correlation container");

    std::vector<AxisSpec> corrAxisTpcFt0c = {{axisSample, "Sample"},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisPtAssoc, "p_{T} (GeV/c)"},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaTpcFt0c, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {
      {axisEtaEfficiency, "#eta"},
      {axisPtEfficiency, "p_{T} (GeV/c)"},
      {axisVertexEfficiency, "z-vtx (cm)"},
    };

    std::vector<AxisSpec> userAxis;

    // Correlation axis For TPC-FIT

    std::vector<AxisSpec> corrAxisTpcFt0a = {{axisSample, "Sample"},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisPtAssoc, "p_{T} (GeV/c)"},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaTpcFt0a, "#Delta#eta"}};
    std::vector<AxisSpec> corrAxisTpcFv0 = {{axisSample, "Sample"},
                                            {axisVertex, "z-vtx (cm)"},
                                            {axisPtTrigger, "p_{T} (GeV/c)"},
                                            {axisPtAssoc, "p_{T} (GeV/c)"},
                                            {axisDeltaPhi, "#Delta#varphi (rad)"},
                                            {axisDeltaEtaTpcFv0, "#Delta#eta"}};

    // correlation axis for TPC-TPC
    std::vector<AxisSpec> corrAxisTpcTpc = {{axisSample, "Sample"},
                                            {axisVertex, "z-vtx (cm)"},
                                            {axisPtTrigger, "p_{T} (GeV/c)"},
                                            {axisPtAssoc, "p_{T} (GeV/c)"},
                                            {axisDeltaPhi, "#Delta#varphi (rad)"},
                                            {axisDeltaEta, "#Delta#eta"}};

    // Correlation axis For TPC-MFT
    std::vector<AxisSpec> corrAxisTpcMft = {{axisSample, "Sample"},
                                            {axisVertex, "z-vtx (cm)"},
                                            {axisPtTrigger, "p_{T} (GeV/c)"},
                                            {axisPtAssoc, "p_{T} (GeV/c)"},
                                            {axisDeltaPhi, "#Delta#varphi (rad)"},
                                            {axisDeltaEtaTpcMft, "#Delta#eta"}};

    // Correlation axis For MFT-FIT
    std::vector<AxisSpec> corrAxisMftFt0a = {{axisSample, "Sample"},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisPtAssoc, "p_{T} (GeV/c)"},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaMftFt0a, "#Delta#eta"}};
    std::vector<AxisSpec> corrAxisMftFt0c = {{axisSample, "Sample"},
                                             {axisVertex, "z-vtx (cm)"},
                                             {axisPtTrigger, "p_{T} (GeV/c)"},
                                             {axisPtAssoc, "p_{T} (GeV/c)"},
                                             {axisDeltaPhi, "#Delta#varphi (rad)"},
                                             {axisDeltaEtaMftFt0c, "#Delta#eta"}};
    std::vector<AxisSpec> corrAxisMftFv0 = {{axisSample, "Sample"},
                                            {axisVertex, "z-vtx (cm)"},
                                            {axisPtTrigger, "p_{T} (GeV/c)"},
                                            {axisPtAssoc, "p_{T} (GeV/c)"},
                                            {axisDeltaPhi, "#Delta#varphi (rad)"},
                                            {axisDeltaEtaMftFv0, "#Delta#eta"}};

    if (doprocessSameTpcFIT) {
      if (cfgDetectorConfig.processFT0A) {
        same.setObject(new CorrelationContainer("sameEvent_TPC_FT0A", "sameEvent_TPC_FT0A", corrAxisTpcFt0a, effAxis, userAxis));
        mixed.setObject(new CorrelationContainer("mixedEvent_TPC_FT0A", "mixedEvent_TPC_FT0A", corrAxisTpcFt0a, effAxis, userAxis));
      }
      if (cfgDetectorConfig.processFT0C) {
        same.setObject(new CorrelationContainer("sameEvent_TPC_FT0C", "sameEvent_TPC_FT0C", corrAxisTpcFt0c, effAxis, userAxis));
        mixed.setObject(new CorrelationContainer("mixedEvent_TPC_FT0C", "mixedEvent_TPC_FT0C", corrAxisTpcFt0c, effAxis, userAxis));
      }
      if (cfgDetectorConfig.processFV0) {
        same.setObject(new CorrelationContainer("sameEvent_TPC_FV0", "sameEvent_TPC_FV0", corrAxisTpcFv0, effAxis, userAxis));
        mixed.setObject(new CorrelationContainer("mixedEvent_TPC_FV0", "mixedEvent_TPC_FV0", corrAxisTpcFv0, effAxis, userAxis));
      }
    }

    if (doprocessSameMFTFIT) {
      if (cfgDetectorConfig.processFT0A) {
        same.setObject(new CorrelationContainer("sameEvent_MFT_FT0A", "sameEvent_MFT_FT0A", corrAxisMftFt0a, effAxis, userAxis));
        mixed.setObject(new CorrelationContainer("mixedEvent_MFT_FT0A", "mixedEvent_MFT_FT0A", corrAxisMftFt0a, effAxis, userAxis));
      }
      if (cfgDetectorConfig.processFT0C) {
        same.setObject(new CorrelationContainer("sameEvent_MFT_FT0C", "sameEvent_MFT_FT0C", corrAxisMftFt0c, effAxis, userAxis));
        mixed.setObject(new CorrelationContainer("mixedEvent_MFT_FT0C", "mixedEvent_MFT_FT0C", corrAxisMftFt0c, effAxis, userAxis));
      }
      if (cfgDetectorConfig.processFV0) {
        same.setObject(new CorrelationContainer("sameEvent_MFT_FV0", "sameEvent_MFT_FV0", corrAxisMftFv0, effAxis, userAxis));
        mixed.setObject(new CorrelationContainer("mixedEvent_MFT_FV0", "mixedEvent_MFT_FV0", corrAxisMftFv0, effAxis, userAxis));
      }
    }

    if (doprocessSameTPC) {
      same.setObject(new CorrelationContainer("sameEvent_TPC_TPC", "sameEvent_TPC_TPC", corrAxisTpcTpc, effAxis, userAxis));
      mixed.setObject(new CorrelationContainer("mixedEvent_TPC_TPC", "mixedEvent_TPC_TPC", corrAxisTpcTpc, effAxis, userAxis));
    }

    if (doprocessSameTpcMft) {
      same.setObject(new CorrelationContainer("sameEvent_TPC_MFT", "sameEvent_TPC_MFT", corrAxisTpcMft, effAxis, userAxis));
      mixed.setObject(new CorrelationContainer("mixedEvent_TPC_MFT", "mixedEvent_TPC_MFT", corrAxisTpcMft, effAxis, userAxis));
    }

    if (doprocessSameTpcMftReassociated2D) {
      same.setObject(new CorrelationContainer("sameEvent_TPC_MFT_Reassociated2D", "sameEvent_TPC_MFT_Reassociated2D", corrAxisTpcMft, effAxis, userAxis));
      mixed.setObject(new CorrelationContainer("mixedEvent_TPC_MFT_Reassociated2D", "mixedEvent_TPC_MFT_Reassociated2D", corrAxisTpcMft, effAxis, userAxis));
    }

    if (doprocessSameTpcMftReassociated3D) {
      same.setObject(new CorrelationContainer("sameEvent_TPC_MFT_Reassociated3D", "sameEvent_TPC_MFT_Reassociated3D", corrAxisTpcMft, effAxis, userAxis));
      mixed.setObject(new CorrelationContainer("mixedEvent_TPC_MFT_Reassociated3D", "mixedEvent_TPC_MFT_Reassociated3D", corrAxisTpcMft, effAxis, userAxis));
    }
    LOGF(info, "End of init");
  }

  TRandom3* gRandom = new TRandom3();

  template <typename TCollision>
  bool eventSelected(TCollision collision, const int multTrk, const bool fillCounter)
  {
    registry.fill(HIST("hEventCountSpecific"), 0.5);
    if (cfgEventSelection.cfgEvSelkNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      // rejects collisions which are associated with the same "found-by-T0" bunch crossing
      // https://indico.cern.ch/event/1396220/#1-event-selection-with-its-rof
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkNoSameBunchPileup)
      registry.fill(HIST("hEventCountSpecific"), 1.5);
    if (cfgEventSelection.cfgEvSelkNoITSROFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkNoITSROFrameBorder)
      registry.fill(HIST("hEventCountSpecific"), 2.5);
    if (cfgEventSelection.cfgEvSelkNoTimeFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkNoTimeFrameBorder)
      registry.fill(HIST("hEventCountSpecific"), 3.5);
    if (cfgEventSelection.cfgEvSelkIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      // removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference
      // use this cut at low multiplicities with caution
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkIsGoodZvtxFT0vsPV)
      registry.fill(HIST("hEventCountSpecific"), 4.5);
    if (cfgEventSelection.cfgEvSelkNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      // no collisions in specified time range
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkNoCollInTimeRangeStandard)
      registry.fill(HIST("hEventCountSpecific"), 5.5);

    if (cfgEventSelection.cfgEvSelkIsGoodITSLayer0123 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayer0123)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkIsGoodITSLayer0123)
      registry.fill(HIST("hEventCountSpecific"), 6.5);

    if (cfgEventSelection.cfgEvSelkIsGoodITSLayersAll && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      // from Jan 9 2025 AOT meeting
      // cut time intervals with dead ITS staves
      return 0;
    }

    if (fillCounter && cfgEventSelection.cfgEvSelkIsGoodITSLayersAll)
      registry.fill(HIST("hEventCountSpecific"), 7.5);

    if (cfgEventSelection.cfgEvSelkNoCollInRofStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      // no other collisions in this Readout Frame with per-collision multiplicity above threshold
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkNoCollInRofStandard)
      registry.fill(HIST("hEventCountSpecific"), 8.5);
    if (cfgEventSelection.cfgEvSelkNoHighMultCollInPrevRof && !collision.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      // veto an event if FT0C amplitude in previous ITS ROF is above threshold
      return 0;
    }
    if (fillCounter && cfgEventSelection.cfgEvSelkNoHighMultCollInPrevRof)
      registry.fill(HIST("hEventCountSpecific"), 9.5);
    auto occupancy = collision.trackOccupancyInTimeRange();
    if (cfgEventSelection.cfgEvSelOccupancy && (occupancy < cfgEventSelection.cfgCutOccupancyLow || occupancy > cfgEventSelection.cfgCutOccupancyHigh))
      return 0;
    if (fillCounter && cfgEventSelection.cfgEvSelOccupancy)
      registry.fill(HIST("hEventCountSpecific"), 10.5);

    auto multNTracksPV = collision.multNTracksPV();

    if (cfgFuncParas.cfgMultGlobalPVCutEnabled) {
      if (multTrk < cfgFuncParas.fMultGlobalPVCutLow->Eval(multNTracksPV))
        return 0;
      if (multTrk > cfgFuncParas.fMultGlobalPVCutHigh->Eval(multNTracksPV))
        return 0;
    }
    if (cfgFuncParas.cfgMultMultV0ACutEnabled) {
      if (collision.multFV0A() < cfgFuncParas.fMultMultV0ACutLow->Eval(multTrk))
        return 0;
      if (collision.multFV0A() > cfgFuncParas.fMultMultV0ACutHigh->Eval(multTrk))
        return 0;
    }

    if (fillCounter && cfgEventSelection.cfgEvSelMultCorrelation)
      registry.fill(HIST("hEventCountSpecific"), 11.5);

    // V0A T0A 5 sigma cut
    float sigma = 5.0;
    if (cfgEventSelection.cfgEvSelV0AT0ACut && (std::fabs(collision.multFV0A() - cfgFuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > sigma * cfgFuncParas.fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;
    if (fillCounter && cfgEventSelection.cfgEvSelV0AT0ACut)
      registry.fill(HIST("hEventCountSpecific"), 12.5);

    return 1;
  }

  double getPhiFV0(uint64_t chno)
  {
    o2::fv0::Point3Dsimple chPos{};
    int const cellsInLeft[] = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 32, 40, 33, 41, 34, 42, 35, 43};
    bool const isChnoInLeft = std::find(std::begin(cellsInLeft), std::end(cellsInLeft), chno) != std::end(cellsInLeft);

    if (isChnoInLeft) {
      chPos = fv0Det->getReadoutCenter(chno);
      return RecoDecay::phi(chPos.x + (*offsetFV0)[0].getX(), chPos.y + (*offsetFV0)[0].getY());
    } else {
      chPos = fv0Det->getReadoutCenter(chno);
      return RecoDecay::phi(chPos.x + (*offsetFV0)[1].getX(), chPos.y + (*offsetFV0)[1].getY());
    }
  }

  double getPhiFT0(uint64_t chno, int i)
  {
    // offsetFT0[0]: FT0A, offsetFT0[1]: FT0C
    if (i > 1 || i < 0) {
      LOGF(fatal, "kFIT Index %d out of range", i);
    }

    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    return RecoDecay::phi(chPos.X() + (*offsetFT0)[i].getX(), chPos.Y() + (*offsetFT0)[i].getY());
  }

  double getEtaFV0(uint64_t chno)
  {

    int const cellsInLeft[] = {0, 1, 2, 3, 8, 9, 10, 11, 16, 17, 18, 19, 24, 25, 26, 27, 32, 40, 33, 41, 34, 42, 35, 43};
    bool const isChnoInLeft = std::find(std::begin(cellsInLeft), std::end(cellsInLeft), chno) != std::end(cellsInLeft);

    o2::fv0::Point3Dsimple chPos{};
    chPos = fv0Det->getReadoutCenter(chno);

    float offsetX, offsetY, offsetZ;

    if (isChnoInLeft) {
      offsetX = (*offsetFV0)[0].getX();
      offsetY = (*offsetFV0)[0].getY();
      offsetZ = (*offsetFV0)[0].getZ();
    } else {
      offsetX = (*offsetFV0)[1].getX();
      offsetY = (*offsetFV0)[1].getY();
      offsetZ = (*offsetFV0)[1].getZ();
    }

    auto x = chPos.x + offsetX;
    auto y = chPos.y + offsetY;
    auto z = chPos.z + offsetZ;

    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);

    return -std::log(std::tan(0.5 * theta));
  }

  double getEtaFT0(uint64_t chno, int i)
  {
    // offsetFT0[0]: FT0A, offsetFT0[1]: FT0C
    if (i > 1 || i < 0) {
      LOGF(fatal, "kFIT Index %d out of range", i);
    }
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    auto x = chPos.X() + (*offsetFT0)[i].getX();
    auto y = chPos.Y() + (*offsetFT0)[i].getY();
    auto z = chPos.Z() + (*offsetFT0)[i].getZ();
    if (chno >= Ft0IndexA) {
      z = -z;
    }
    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);
    return -std::log(std::tan(0.5 * theta));
  }

  template <typename TTrack>
  bool isAcceptedMftTrack(TTrack const& mftTrack)
  {
    // cut on the eta of MFT tracks
    if (mftTrack.eta() < cfgMftConfig.etaMftTrackMin || mftTrack.eta() > cfgMftConfig.etaMftTrackMax) {
      return false;
    }

    // cut on the number of clusters of the reconstructed MFT track
    if (mftTrack.nClusters() < cfgMftConfig.nClustersMftTrack) {
      return false;
    }

    if (mftTrack.pt() < cfgMftConfig.cfgPtCutMinMFT || mftTrack.pt() > cfgMftConfig.cfgPtCutMaxMFT) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isAmbiguousMftTrack(TTrack const& mftTrack, bool fillHistogram)
  {
    if (mftTrack.ambDegree() > 1) {
      if (fillHistogram) {
        registry.fill(HIST("hEventCountMftReassoc"), 2.5); // fill histogram for events with at least one ambiguous track);
      }
      return false;
    }
    registry.fill(HIST("hEventCountMftReassoc"), 3.5); // fill histogram for events without ambiguous tracks
    return true;
  }

  void loadAlignParam(uint64_t timestamp)
  {
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", timestamp);
    offsetFV0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FV0/Calib/Align", timestamp);
    if (offsetFT0 == nullptr) {
      LOGF(fatal, "Could not load FT0/Calib/Align for timestamp %d", timestamp);
    }
    if (offsetFV0 == nullptr) {
      LOGF(fatal, "Could not load FV0/Calib/Align for timestamp %d", timestamp);
    }
  }

  void loadGain(aod::BCsWithTimestamps::iterator const& bc)
  {
    cstFT0RelGain.clear();
    cstFT0RelGain = {};
    std::string fullPath;

    auto timestamp = bc.timestamp();
    constexpr int ChannelsFT0 = 208;
    if (cfgCorrLevel == 0) {
      for (auto i{0u}; i < ChannelsFT0; i++) {
        cstFT0RelGain.push_back(1.);
      }
    } else {
      fullPath = cfgGainEqPath;
      fullPath += "/FT0";
      const auto objft0Gain = ccdb->getForTimeStamp<std::vector<float>>(fullPath, timestamp);
      if (!objft0Gain) {
        for (auto i{0u}; i < ChannelsFT0; i++) {
          cstFT0RelGain.push_back(1.);
        }
      } else {
        cstFT0RelGain = *(objft0Gain);
      }
    }
  }

  void loadCorrection(uint64_t timestamp)
  {
    if (correctionsLoaded) {
      return;
    }
    if (cfgEfficiency.value.empty() == false) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyTrigger = TFile::Open(cfgEfficiency.value.c_str(), "READ");
        mEfficiency = reinterpret_cast<TH3D*>(fEfficiencyTrigger->Get("ccdb_object"));
      } else {
        mEfficiency = ccdb->getForTimeStamp<TH3D>(cfgEfficiency, timestamp);
      }
      if (mEfficiency == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiency.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgEfficiency.value.c_str(), (void*)mEfficiency);
    }
    if (cfgCentralityWeight.value.empty() == false) {
      mCentralityWeight = ccdb->getForTimeStamp<TH1D>(cfgCentralityWeight, timestamp);
      if (mCentralityWeight == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgCentralityWeight.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram from %s (%p)", cfgCentralityWeight.value.c_str(), (void*)mCentralityWeight);
    }
    correctionsLoaded = true;
  }

  bool getEfficiencyCorrection(float& weight_nue, float eta, float pt, float posZ)
  {
    float eff = 1.;
    if (mEfficiency) {
      int etaBin = mEfficiency->GetXaxis()->FindBin(eta);
      int ptBin = mEfficiency->GetYaxis()->FindBin(pt);
      int zBin = mEfficiency->GetZaxis()->FindBin(posZ);
      eff = mEfficiency->GetBinContent(etaBin, ptBin, zBin);
    } else {
      eff = 1.0;
    }
    if (eff == 0)
      return false;
    weight_nue = 1. / eff;
    return true;
  }

  template <typename TFT0s>
  void getChannelFT0(TFT0s const& ft0, std::size_t const& iCh, int& id, float& ampl, int fitType)
  {
    if (fitType == kFT0C) {
      id = ft0.channelC()[iCh];
      id = id + Ft0IndexA;
      ampl = ft0.amplitudeC()[iCh];
    } else if (fitType == kFT0A) {
      id = ft0.channelA()[iCh];
      ampl = ft0.amplitudeA()[iCh];
    } else {
      LOGF(fatal, "Cor Index %d out of range", fitType);
    }
    registry.fill(HIST("FT0Amp"), id, ampl);
  }

  template <typename TFT0s>
  void getChannelFV0(TFT0s const& fv0, std::size_t const& iCh, int& id, float& ampl)
  {
    id = fv0.channel()[iCh];
    ampl = fv0.amplitude()[iCh];
    registry.fill(HIST("FV0Amp"), id, ampl);
  }

  template <typename TTrack>
  bool trackSelected(TTrack track)
  {
    return ((track.tpcNClsFound() >= cfgTrackCuts.cfgCutTPCclu) && (track.tpcNClsCrossedRows() >= cfgTrackCuts.cfgCutTPCCrossedRows) && (track.itsNCls() >= cfgTrackCuts.cfgCutITSclu));
  }

  int getMagneticField(uint64_t timestamp)
  {
    // Get the magnetic field
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("/GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  // fill multiple histograms
  template <typename TCollision, typename TTracks>
  void fillYield(TCollision collision, TTracks tracks) // function to fill the yield and etaphi histograms.
  {
    registry.fill(HIST("zVtx"), collision.posZ());
    registry.fill(HIST("Nch"), tracks.size());

    float weff1 = 1.0;
    float zvtx = collision.posZ();

    for (auto const& track1 : tracks) {

      if constexpr (std::is_same_v<aod::MFTTracks, TTracks>) {
        if (!isAcceptedMftTrack(track1)) {
          continue;
        }
      } else {
        if (!trackSelected(track1)) {
          continue;
        }
        if (!getEfficiencyCorrection(weff1, track1.eta(), track1.pt(), zvtx)) {
          continue;
        }
      }

      registry.fill(HIST("Phi"), RecoDecay::constrainAngle(track1.phi(), 0.0));
      registry.fill(HIST("Eta"), track1.eta());
      registry.fill(HIST("EtaCorrected"), track1.eta(), weff1);
      registry.fill(HIST("pT"), track1.pt());
      registry.fill(HIST("pTCorrected"), track1.pt(), weff1);
    }
  }

  template <typename TTrack, typename TTrackAssoc>
  float getDPhiStar(TTrack const& track1, TTrackAssoc const& track2, float radius, int magField)
  {
    float charge1 = track1.sign();
    float charge2 = track2.sign();

    float phi1 = track1.phi();
    float phi2 = track2.phi();

    float pt1 = track1.pt();
    float pt2 = track2.pt();

    int fbSign = (magField > 0) ? 1 : -1;

    float dPhiStar = phi1 - phi2 - charge1 * fbSign * std::asin(0.075 * radius / pt1) + charge2 * fbSign * std::asin(0.075 * radius / pt2);

    if (dPhiStar > constants::math::PI)
      dPhiStar = constants::math::TwoPI - dPhiStar;
    if (dPhiStar < -constants::math::PI)
      dPhiStar = -constants::math::TwoPI - dPhiStar;

    return dPhiStar;
  }

  // Correlations for detectors and TPC

  //////////////////////////////
  //////////MFT/TPC-FIT////////
  ////////////////////////////
  template <CorrelationContainer::CFStep step, typename TTracks, typename TTracksAssociated, typename FITs>
  void fillCorrelationsFIT(TTracks tracks1, TTracksAssociated tracks2, FITs const&, float posZ, int system, int corType, float multiplicity)
  {

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    float triggerWeight = 1.0f;

    // loop over all tracks

    if (system == SameEvent) {
      registry.fill(HIST("Nch_used"), multiplicity);
    }

    for (auto const& track1 : tracks1) {

      if constexpr (std::is_same_v<aod::MFTTracks, TTracks>) {

        if (!isAcceptedMftTrack(track1)) {
          continue;
        }
      } else {
        if (!trackSelected(track1)) {
          continue;
        }

        if (!getEfficiencyCorrection(triggerWeight, track1.eta(), track1.pt(), posZ))
          continue;
      }

      if (system == SameEvent) {
        registry.fill(HIST("Trig_hist"), fSampleIndex, posZ, track1.pt(), triggerWeight);
      }

      // if using FV0A for correlations / using FV0A as associated particles
      if constexpr (std::is_same_v<aod::FV0As, FITs>) {

        std::size_t channelSize = tracks2.channel().size();
        for (std::size_t iCh = 0; iCh < channelSize; iCh++) {
          int channelID = 0;
          float amplitude = 0.;

          getChannelFV0(tracks2, iCh, channelID, amplitude);

          auto phi = getPhiFV0(channelID);
          auto eta = getEtaFV0(channelID);

          float deltaPhi = RecoDecay::constrainAngle(track1.phi() - phi, -PIHalf);
          float deltaEta = track1.eta() - eta;

          if (system == SameEvent) {
            registry.fill(HIST("Assoc_amp_same"), channelID, amplitude);
            same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track1.pt(), deltaPhi, deltaEta, amplitude * triggerWeight);
            registry.fill(HIST("deltaEta_deltaPhi_same_MFT_FV0"), deltaPhi, deltaEta, amplitude * triggerWeight);
          } else if (system == MixedEvent) {
            registry.fill(HIST("Assoc_amp_mixed"), channelID, amplitude);
            mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track1.pt(), deltaPhi, deltaEta, amplitude);
            registry.fill(HIST("deltaEta_deltaPhi_mixed_MFT_FV0"), deltaPhi, deltaEta, amplitude);
          }
        }
      }

      // if using FT0A and FT0C for correlations / using FT0A and FT0C as associated particles
      if constexpr (std::is_same_v<aod::FT0s, FITs>) {

        std::size_t channelSize = 0;
        if (corType == kFT0C) {
          channelSize = tracks2.channelC().size();
        } else if (corType == kFT0A) {
          channelSize = tracks2.channelA().size();
        } else {
          LOGF(fatal, "Cor Index %d out of range", corType);
        }

        for (std::size_t iCh = 0; iCh < channelSize; iCh++) {
          int channelID = 0;
          float amplitude = 0.;
          getChannelFT0(tracks2, iCh, channelID, amplitude, corType);

          // reject depending on FT0C/FT0A rings
          if (corType == kFT0C) {
            if ((cfgFITConfig.cfgRejectFT0CInside && (channelID >= kFT0CInnerRingMin && channelID <= kFT0CInnerRingMax)) || (cfgFITConfig.cfgRejectFT0COutside && (channelID >= kFT0COuterRingMin && channelID <= kFT0COuterRingMax)))
              continue;
          }
          if (corType == kFT0A) {
            if ((cfgFITConfig.cfgRejectFT0AInside && (channelID >= kFT0AInnerRingMin && channelID <= kFT0AInnerRingMax)) || (cfgFITConfig.cfgRejectFT0AOutside && (channelID >= kFT0AOuterRingMin && channelID <= kFT0AOuterRingMax)))
              continue;
          }

          auto phi = getPhiFT0(channelID, corType);
          auto eta = getEtaFT0(channelID, corType);

          float deltaPhi = RecoDecay::constrainAngle(track1.phi() - phi, -PIHalf);
          float deltaEta = track1.eta() - eta;

          if (system == SameEvent) {
            if (corType == kFT0A) {
              registry.fill(HIST("Assoc_amp_same"), channelID, amplitude);
              same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track1.pt(), deltaPhi, deltaEta, amplitude * triggerWeight);
              registry.fill(HIST("deltaEta_deltaPhi_same_MFT_FT0A"), deltaPhi, deltaEta, amplitude * triggerWeight);
            }
            if (corType == kFT0C) {
              registry.fill(HIST("Assoc_amp_same"), channelID, amplitude);
              same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track1.pt(), deltaPhi, deltaEta, amplitude * triggerWeight);
              registry.fill(HIST("deltaEta_deltaPhi_same_MFT_FT0C"), deltaPhi, deltaEta, amplitude * triggerWeight);
            }
          } else if (system == MixedEvent) {
            if (corType == kFT0A) {
              registry.fill(HIST("Assoc_amp_mixed"), channelID, amplitude);
              mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track1.pt(), deltaPhi, deltaEta, amplitude);
              registry.fill(HIST("deltaEta_deltaPhi_mixed_MFT_FT0A"), deltaPhi, deltaEta, amplitude);
            }
            if (corType == kFT0C) {
              registry.fill(HIST("Assoc_amp_mixed"), channelID, amplitude);
              mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track1.pt(), deltaPhi, deltaEta, amplitude);
              registry.fill(HIST("deltaEta_deltaPhi_mixed_MFT_FT0C"), deltaPhi, deltaEta, amplitude);
            }
          }
        }
      }
    }
  }
  //////////////////////////
  //////////TPC-MFT/////////
  //////////////////////////
  template <CorrelationContainer::CFStep step, typename TTracks, typename TTracksAssoc>
  void fillCorrelationsMFT(TTracks tracks1, TTracksAssoc tracks2, float posZ, int system, int magneticField) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    float triggerWeight = 1.0f;

    if (system == SameEvent) {
      registry.fill(HIST("Nch_used"), tracks1.size());
    }
    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1))
        continue;

      if (!getEfficiencyCorrection(triggerWeight, track1.eta(), track1.pt(), posZ))
        continue;

      if (system == SameEvent) {
        registry.fill(HIST("Trig_hist_TPC_MFT"), fSampleIndex, posZ, track1.pt(), triggerWeight);
      }

      for (auto const& track2 : tracks2) {
        if constexpr (std::is_same_v<aod::MFTTracks, TTracksAssoc>)
          continue;

        if (!isAcceptedMftTrack(track2)) {
          continue;
        }

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        if (cfgApplyTwoTrackEfficiency && std::abs(deltaEta) < cfgMergingCut) {

          double dPhiStarHigh = getDPhiStar(track1, track2, cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(track1, track2, cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgMergingCut;

          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgRadiusLow); rad < cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getDPhiStar(track1, track2, rad, magneticField);
              if (std::abs(dPhiStar) < kLimit) {
                bIsBelow = true;
                break;
              }
            }
            if (bIsBelow)
              continue;
          }
        }

        // fill the right sparse and histograms
        if (system == SameEvent) {

          same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta);
          registry.fill(HIST("deltaEta_deltaPhi_same_TPC_MFT"), deltaPhi, deltaEta);
        } else if (system == MixedEvent) {

          mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta);
          registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC_MFT"), deltaPhi, deltaEta);
        }
      }
    }
  }

  template <CorrelationContainer::CFStep step, typename TTracks, typename TTracksAssoc>
  void fillCorrelationsMftReassociatedTracks(TTracks tracks1, TTracksAssoc tracks2, float multiplicity, float posZ, int system, int magneticField, bool cutAmbiguousTracks) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    float triggerWeight = 1.0f;

    auto loopCounter = 0;

    if (system == SameEvent) {
      registry.fill(HIST("Nch_used"), multiplicity);
    }

    // loop over all tracks
    for (auto const& track1 : tracks1) {

      loopCounter++;

      if (!trackSelected(track1))
        continue;

      if (!getEfficiencyCorrection(triggerWeight, track1.eta(), track1.pt(), posZ))
        continue;

      if (system == SameEvent) {
        registry.fill(HIST("Trig_hist_TPC_MFT"), fSampleIndex, posZ, track1.pt(), triggerWeight);
      }

      for (auto const& track2 : tracks2) {

        auto reassociatedMftTrack = track2.template mfttrack_as<FilteredMftTracks>();

        if (!cutAmbiguousTracks && system == SameEvent && (loopCounter == 1)) {
          registry.fill(HIST("hEventCountMftReassoc"), 0.5); // fill histogram for events with at least one reassociated track);
        }

        if (!isAcceptedMftTrack(reassociatedMftTrack)) {
          continue;
        }

        if (!cutAmbiguousTracks && system == SameEvent && (loopCounter == 1)) {
          registry.fill(HIST("hEventCountMftReassoc"), 1.5); // fill histogram for events with at least one reassociated track after track selection);
        }

        if (isAmbiguousMftTrack(track2, (!cutAmbiguousTracks && system == SameEvent && (loopCounter == 1)))) {

          if (SameEvent && (loopCounter == 1)) {
            registry.fill(HIST("ReassociatedMftTracks"), 0.5);
          }
          if (cutAmbiguousTracks) {
            continue;
          }
        }

        if (reassociatedMftTrack.collisionId() != track2.bestCollisionId()) {
          if (SameEvent && (loopCounter == 1)) {
            registry.fill(HIST("ReassociatedMftTracks"), 1.5);
          }
        }

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - reassociatedMftTrack.phi(), -PIHalf);
        float deltaEta = track1.eta() - reassociatedMftTrack.eta();

        if (cfgApplyTwoTrackEfficiency && std::abs(deltaEta) < cfgMergingCut) {

          double dPhiStarHigh = getDPhiStar(track1, reassociatedMftTrack, cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(track1, reassociatedMftTrack, cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgMergingCut;

          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgRadiusLow); rad < cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getDPhiStar(track1, reassociatedMftTrack, rad, magneticField);
              if (std::abs(dPhiStar) < kLimit) {
                bIsBelow = true;
                break;
              }
            }
            if (bIsBelow)
              continue;
          }
        }

        // fill the right sparse and histograms
        if (system == SameEvent) {

          same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), reassociatedMftTrack.pt(), deltaPhi, deltaEta);
          registry.fill(HIST("deltaEta_deltaPhi_same_TPC_MFT"), deltaPhi, deltaEta);
        } else if (system == MixedEvent) {

          mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), reassociatedMftTrack.pt(), deltaPhi, deltaEta);
          registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC_MFT"), deltaPhi, deltaEta);
        }
      }
    }
  }

  ///////////////////////////////////////
  //////////TPC-TPC and TPC-MFT/////////
  /////////////////////////////////////
  template <CorrelationContainer::CFStep step, typename TTracks, typename TTracksAssoc>
  void fillCorrelationsTpc(TTracks tracks1, TTracksAssoc tracks2, float posZ, int system, int magneticField) // function to fill the Output functions (sparse) and the delta eta and delta phi histograms
  {

    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);

    float triggerWeight = 1.0f;

    // loop over all tracks
    for (auto const& track1 : tracks1) {

      if (!trackSelected(track1))
        continue;

      if (!getEfficiencyCorrection(triggerWeight, track1.eta(), track1.pt(), posZ))
        continue;

      if (system == SameEvent) {
        registry.fill(HIST("Nch_used"), tracks1.size());
        registry.fill(HIST("Trig_hist_TPC_TPC"), fSampleIndex, posZ, track1.pt(), triggerWeight);
      }

      for (auto const& track2 : tracks2) {

        if (cfgDetectorConfig.processMFT) {
          if constexpr (std::is_same_v<aod::MFTTracks, TTracksAssoc>) {
            if (!isAcceptedMftTrack(track2)) {
              continue;
            }
          }
        } else {
          if (!trackSelected(track2))
            continue;

          if (track1.pt() <= track2.pt())
            continue; // skip if the trigger pt is less than the associate pt
        }

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -PIHalf);
        float deltaEta = track1.eta() - track2.eta();

        if (cfgApplyTwoTrackEfficiency && std::abs(deltaEta) < cfgMergingCut) {

          double dPhiStarHigh = getDPhiStar(track1, track2, cfgRadiusHigh, magneticField);
          double dPhiStarLow = getDPhiStar(track1, track2, cfgRadiusLow, magneticField);

          const double kLimit = 3.0 * cfgMergingCut;

          bool bIsBelow = false;

          if (std::abs(dPhiStarLow) < kLimit || std::abs(dPhiStarHigh) < kLimit || dPhiStarLow * dPhiStarHigh < 0) {
            for (double rad(cfgRadiusLow); rad < cfgRadiusHigh; rad += 0.01) {
              double dPhiStar = getDPhiStar(track1, track2, rad, magneticField);
              if (std::abs(dPhiStar) < kLimit) {
                bIsBelow = true;
                break;
              }
            }
            if (bIsBelow)
              continue;
          }
        }

        // fill the right sparse and histograms
        if (system == SameEvent) {
          if (cfgDetectorConfig.processMFT) {
            same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta);
            registry.fill(HIST("deltaEta_deltaPhi_same_TPC_MFT"), deltaPhi, deltaEta);
          } else {
            same->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta);
            registry.fill(HIST("deltaEta_deltaPhi_same_TPC_TPC"), deltaPhi, deltaEta);
          }
        } else if (system == MixedEvent) {
          if (cfgDetectorConfig.processMFT) {
            mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta);
            registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC_MFT"), deltaPhi, deltaEta);
          } else {
            mixed->getPairHist()->Fill(step, fSampleIndex, posZ, track1.pt(), track2.pt(), deltaPhi, deltaEta);
            registry.fill(HIST("deltaEta_deltaPhi_mixed_TPC_TPC"), deltaPhi, deltaEta);
          }
        }
      }
    }
  }

  //////////////////////////////////////
  ////// Same event processing /////////
  //////////////////////////////////////

  ////////////////////////////////////
  /////////// Fwrd-Bwrd /////////////
  ////////////////////////////////////

  void processSameMFTFIT(AodCollisions::iterator const& collision, AodTracks const& tpctracks, aod::MFTTracks const& mfts, aod::FT0s const& ft0as, aod::FV0As const& fv0as, aod::BCsWithTimestamps const&)
  {

    if (!collision.sel8())
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (cfgUseAdditionalEventCut && !eventSelected(collision, tpctracks.size(), true))
      return;

    if (!collision.has_foundFT0())
      return;
    loadAlignParam(bc.timestamp());
    // loadGain(bc);
    loadCorrection(bc.timestamp());

    if ((tpctracks.size() < cfgEventSelection.cfgMinMult || tpctracks.size() >= cfgEventSelection.cfgMaxMult)) {
      return;
    }

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin
    fillYield(collision, mfts);
    const auto& multiplicity = tpctracks.size();

    if (cfgDetectorConfig.processFV0) {
      if (collision.has_foundFV0()) {
        same->fillEvent(mfts.size(), CorrelationContainer::kCFStepReconstructed);
        const auto& fv0 = collision.foundFV0();
        fillCorrelationsFIT<CorrelationContainer::kCFStepReconstructed>(mfts, fv0, fv0as, collision.posZ(), SameEvent, kFV0, multiplicity);
      }
    }
    if (cfgDetectorConfig.processFT0C) {
      if (collision.has_foundFT0()) {
        same->fillEvent(mfts.size(), CorrelationContainer::kCFStepReconstructed);
        const auto& ft0 = collision.foundFT0();
        fillCorrelationsFIT<CorrelationContainer::kCFStepReconstructed>(mfts, ft0, ft0as, collision.posZ(), SameEvent, kFT0C, multiplicity);
      }
    }
    if (cfgDetectorConfig.processFT0A) {
      if (collision.has_foundFT0()) {
        same->fillEvent(mfts.size(), CorrelationContainer::kCFStepReconstructed);
        const auto& ft0 = collision.foundFT0();
        fillCorrelationsFIT<CorrelationContainer::kCFStepReconstructed>(mfts, ft0, ft0as, collision.posZ(), SameEvent, kFT0A, multiplicity);
      }
    }
  }
  PROCESS_SWITCH(CorrSparse, processSameMFTFIT, "Process same event for MFT-FIT correlation", true);

  /////////////////////////
  ////////Mid-Mid//////////
  ////////////////////////

  void processSameTPC(AodCollisions::iterator const& collision, AodTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (!collision.sel8())
      return;

    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true))
      return;

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    if (tracks.size() < cfgEventSelection.cfgMinMult || tracks.size() >= cfgEventSelection.cfgMaxMult) {
      return;
    }

    loadCorrection(bc.timestamp());
    fillYield(collision, tracks);

    fillCorrelationsTpc<CorrelationContainer::kCFStepReconstructed>(tracks, tracks, collision.posZ(), SameEvent, getMagneticField(bc.timestamp()));
  }
  PROCESS_SWITCH(CorrSparse, processSameTPC, "Process same event for TPC-TPC correlation", false);

  /////////////////////
  ////////back-Mid-Fwrd//////
  /////////////////////

  void processSameTpcFIT(AodCollisions::iterator const& collision, AodTracks const& tracks, aod::FT0s const& ft0as, aod::FV0As const& fv0as, aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8())
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true))
      return;
    if (!collision.has_foundFT0() && !collision.has_foundFV0())
      return;

    registry.fill(HIST("Nch"), tracks.size());
    registry.fill(HIST("zVtx"), collision.posZ());

    if ((tracks.size() < cfgEventSelection.cfgMinMult || tracks.size() >= cfgEventSelection.cfgMaxMult)) {
      return;
    }

    loadAlignParam(bc.timestamp());
    // loadGain(bc);
    loadCorrection(bc.timestamp());

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin
    fillYield(collision, tracks);

    const auto& multiplicity = tracks.size();

    if (cfgDetectorConfig.processFV0) {
      same->fillEvent(tracks.size(), CorrelationContainer::kCFStepReconstructed);
      const auto& fv0 = collision.foundFV0();
      fillCorrelationsFIT<CorrelationContainer::kCFStepReconstructed>(tracks, fv0, fv0as, collision.posZ(), SameEvent, kFV0, multiplicity);
    }
    if (cfgDetectorConfig.processFT0C) {
      same->fillEvent(tracks.size(), CorrelationContainer::kCFStepReconstructed);
      const auto& ft0 = collision.foundFT0();
      fillCorrelationsFIT<CorrelationContainer::kCFStepReconstructed>(tracks, ft0, ft0as, collision.posZ(), SameEvent, kFT0C, multiplicity);
    }
    if (cfgDetectorConfig.processFT0A) {
      same->fillEvent(tracks.size(), CorrelationContainer::kCFStepReconstructed);
      const auto& ft0 = collision.foundFT0();
      fillCorrelationsFIT<CorrelationContainer::kCFStepReconstructed>(tracks, ft0, ft0as, collision.posZ(), SameEvent, kFT0A, multiplicity);
    }
  }
  PROCESS_SWITCH(CorrSparse, processSameTpcFIT, "process for forward or backwards correlations with TPC", false);

  void processSameTpcMft(AodCollisions::iterator const& collision, AodTracks const& tracks, aod::MFTTracks const& mfts, aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8())
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true))
      return;

    loadCorrection(bc.timestamp());
    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    fillYield(collision, tracks);

    if (tracks.size() < cfgEventSelection.cfgMinMult || tracks.size() >= cfgEventSelection.cfgMaxMult) {
      return;
    }

    fillCorrelationsMFT<CorrelationContainer::kCFStepReconstructed>(tracks, mfts, collision.posZ(), SameEvent, getMagneticField(bc.timestamp()));
  }
  PROCESS_SWITCH(CorrSparse, processSameTpcMft, "Process same event for TPC-MFT correlation", false);

  void processSameTpcMftReassociated2D(AodCollisions::iterator const& collision, AodTracks const& tracks,
                                       soa::SmallGroups<aod::BestCollisionsFwd> const& reassociatedMftTracks,
                                       FilteredMftTracks const&,
                                       aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8())
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true))
      return;

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin

    loadCorrection(bc.timestamp());
    fillYield(collision, tracks);

    if (tracks.size() < cfgEventSelection.cfgMinMult || tracks.size() >= cfgEventSelection.cfgMaxMult) {
      return;
    }

    fillCorrelationsMftReassociatedTracks<CorrelationContainer::kCFStepReconstructed>(tracks, reassociatedMftTracks, collision.posZ(), tracks.size(), SameEvent, getMagneticField(bc.timestamp()), true);
  }
  PROCESS_SWITCH(CorrSparse, processSameTpcMftReassociated2D, "Process same event for TPC-MFT correlation with reassociated tracks", false);

  void processSameTpcMftReassociated3D(AodCollisions::iterator const& collision, AodTracks const& tracks,
                                       soa::SmallGroups<aod::BestCollisionsFwd3d> const& reassociatedMftTracks,
                                       aod::BCsWithTimestamps const&)
  {
    if (!collision.sel8())
      return;
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (cfgUseAdditionalEventCut && !eventSelected(collision, tracks.size(), true))
      return;

    registry.fill(HIST("eventcount"), SameEvent); // because its same event i put it in the 1 bin
    loadCorrection(bc.timestamp());
    fillYield(collision, tracks);

    if (tracks.size() < cfgEventSelection.cfgMinMult || tracks.size() >= cfgEventSelection.cfgMaxMult) {
      return;
    }

    fillCorrelationsMftReassociatedTracks<CorrelationContainer::kCFStepReconstructed>(tracks, reassociatedMftTracks, collision.posZ(), tracks.size(), SameEvent, getMagneticField(bc.timestamp()), true);
  }
  PROCESS_SWITCH(CorrSparse, processSameTpcMftReassociated3D, "Process same event for TPC-MFT correlation with reassociated tracks", false);

  ////////////////////////////////////
  ////// Mixed event processing //////
  ////////////////////////////////////

  /////////////////////////////////////////
  ////////////// Fwrd- Bwrd//////////////
  ////////////////////////////////////////

  void processMixedMFTFIT(AodCollisions const& collisions, AodTracks const& tpctracks, aod::MFTTracks const& mfts, aod::FT0s const& ft0as, aod::FV0As const& fv0as, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tpctracks, this](AodCollisions::iterator const& collision) {
      auto associatedTracks = tpctracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {vtxMix, multMix}, true};

    auto tracksTuple = std::make_tuple(mfts, mfts);
    Pair<AodCollisions, aod::MFTTracks, aod::MFTTracks, MixedBinning> pair{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto it = pair.begin(); it != pair.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      auto slicedtracks = tpctracks.sliceBy(perColGlobal, collision1.globalIndex());
      auto multiplicity = slicedtracks.size();

      if ((multiplicity < cfgEventSelection.cfgMinMult || multiplicity >= cfgEventSelection.cfgMaxMult))
        continue;

      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false))
        continue;

      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false))
        continue;

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadAlignParam(bc.timestamp());

      loadCorrection(bc.timestamp());

      if (cfgDetectorConfig.processFT0A) {
        if (!collision1.has_foundFT0() && !collision2.has_foundFT0())
          continue;

        const auto& ft0 = collision2.foundFT0();
        fillCorrelationsFIT<CorrelationContainer::kCFStepReconstructed>(tracks1, ft0, ft0as, collision1.posZ(), MixedEvent, kFT0A, multiplicity);
      }
      if (cfgDetectorConfig.processFT0C) {
        if (!collision1.has_foundFT0() && !collision2.has_foundFT0())
          continue;

        const auto& ft0 = collision2.foundFT0();
        fillCorrelationsFIT<CorrelationContainer::kCFStepReconstructed>(tracks1, ft0, ft0as, collision1.posZ(), MixedEvent, kFT0C, multiplicity);
      }
      if (cfgDetectorConfig.processFV0) {
        if (collision1.has_foundFV0() && collision2.has_foundFV0()) {
          const auto& fv0 = collision2.foundFV0();
          fillCorrelationsFIT<CorrelationContainer::kCFStepReconstructed>(tracks1, fv0, fv0as, collision1.posZ(), MixedEvent, kFV0, multiplicity);
        }
      }
    }
  }
  PROCESS_SWITCH(CorrSparse, processMixedMFTFIT, "Process mixed events for MFT-FIT correlation", true);

  /////////////////////////////
  //////////Mid-Mid///////////
  ////////////////////////////

  void processMixedTpcTpc(AodCollisions const& collisions, AodTracks const& tracks, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](AodCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {vtxMix, multMix}, true};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<AodCollisions, AodTracks, AodTracks, MixedBinning> pair{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto const& [collision1, tracks1, collision2, tracks2] : pair) {

      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();

      loadCorrection(bc.timestamp());

      if ((tracks1.size() < cfgEventSelection.cfgMinMult || tracks1.size() >= cfgEventSelection.cfgMaxMult))
        continue;

      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false))
        continue;
      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false))
        continue;

      fillCorrelationsTpc<CorrelationContainer::kCFStepReconstructed>(tracks1, tracks2, collision1.posZ(), MixedEvent, getMagneticField(bc.timestamp()));
    }
  }
  PROCESS_SWITCH(CorrSparse, processMixedTpcTpc, "Process mixed events for TPC-TPC correlation", false);

  //////////////////////////////
  /////// back-Mid-Fwrd ///////
  /////////////////////////////

  void processMixedTpcFIT(AodCollisions const& collisions, AodTracks const& tracks, aod::FV0As const& fv0as, aod::FT0s const& ft0as, aod::BCsWithTimestamps const&)
  {
    auto getTracksSize = [&tracks, this](AodCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {vtxMix, multMix}, true};

    auto tracksTuple = std::make_tuple(tracks, tracks);
    Pair<AodCollisions, AodTracks, AodTracks, MixedBinning> pair{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto it = pair.begin(); it != pair.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      if ((tracks1.size() < cfgEventSelection.cfgMinMult || tracks1.size() >= cfgEventSelection.cfgMaxMult))
        continue;

      if (cfgUseAdditionalEventCut && !eventSelected(collision1, tracks1.size(), false))
        continue;
      if (cfgUseAdditionalEventCut && !eventSelected(collision2, tracks2.size(), false))
        continue;

      if (!(collision1.has_foundFT0() && collision2.has_foundFT0()))
        continue;

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
      loadAlignParam(bc.timestamp());
      loadCorrection(bc.timestamp());

      auto multiplicity = tracks1.size();

      if (cfgDetectorConfig.processFT0A) {
        const auto& ft0 = collision2.foundFT0();
        fillCorrelationsFIT<CorrelationContainer::kCFStepReconstructed>(tracks1, ft0, ft0as, collision1.posZ(), MixedEvent, kFT0A, multiplicity);
      }
      if (cfgDetectorConfig.processFT0C) {
        const auto& ft0 = collision2.foundFT0();
        fillCorrelationsFIT<CorrelationContainer::kCFStepReconstructed>(tracks1, ft0, ft0as, collision1.posZ(), MixedEvent, kFT0C, multiplicity);
      }
      if (cfgDetectorConfig.processFV0) {
        const auto& fv0 = collision2.foundFV0();
        fillCorrelationsFIT<CorrelationContainer::kCFStepReconstructed>(tracks1, fv0, fv0as, collision1.posZ(), MixedEvent, kFV0, multiplicity);
      }
    }
  }
  PROCESS_SWITCH(CorrSparse, processMixedTpcFIT, "Process mixed events for TPC-FIT correlation", false);

  void processMixedTpcMFT(AodCollisions const& collisions, AodTracks const& tracks, aod::MFTTracks const& MFTtracks, aod::BCsWithTimestamps const&)
  {

    auto getTracksSize = [&tracks, this](AodCollisions::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      auto mult = associatedTracks.size();
      return mult;
    };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;

    MixedBinning binningOnVtxAndMult{{getTracksSize}, {vtxMix, multMix}, true};

    auto tracksTuple = std::make_tuple(tracks, MFTtracks);
    Pair<AodCollisions, AodTracks, aod::MFTTracks, MixedBinning> pair{binningOnVtxAndMult, cfgMinMixEventNum, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip
    for (auto const& [collision1, tracks1, collision2, tracks2] : pair) {

      if (!collision1.sel8() || !collision2.sel8())
        continue;

      if (collision1.globalIndex() == collision2.globalIndex()) {
        continue;
      }

      registry.fill(HIST("eventcount"), MixedEvent); // fill the mixed event in the 3 bin
      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();

      loadCorrection(bc.timestamp());

      if ((tracks1.size() < cfgEventSelection.cfgMinMult || tracks1.size() >= cfgEventSelection.cfgMaxMult))
        continue;

      fillCorrelationsMFT<CorrelationContainer::kCFStepReconstructed>(tracks1, tracks2, collision1.posZ(), MixedEvent, getMagneticField(bc.timestamp()));
    }
  }
  PROCESS_SWITCH(CorrSparse, processMixedTpcMFT, "Process mixed events for TPC-MFT correlation", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CorrSparse>(cfgc),
  };
}
