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

/// \file taskFlowCharmHadrons.cxx
/// \brief Analysis task for charm hadron flow
///
/// \author S. Politanò, INFN Torino, Italy
/// \author Wu Chuntai, CUG, China
/// \author Ran Tu, Fudan University, China
/// \author Marcello Di Costanzo <marcello.di.costanzo@cern.ch>, Polytechnic University of Turin and INFN

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/D2H/Utils/utilsFlow.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>
#include <TString.h>

#include <Rtypes.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;
using namespace o2::hf_evsel;
using namespace o2::analysis::hf_flow_utils;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(M, m, float);                   //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                 //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(MlScore0, mlScore0, float);     //! ML score of the first configured index
DECLARE_SOA_COLUMN(MlScore1, mlScore1, float);     //! ML score of the second configured index
DECLARE_SOA_COLUMN(ScalarProd, scalarProd, float); //! Scalar product
DECLARE_SOA_COLUMN(Cent, cent, float);             //! Centrality
} // namespace full
DECLARE_SOA_TABLE(HfCandMPtInfos, "AOD", "HFCANDMPTINFO",
                  full::M,
                  full::Pt,
                  full::MlScore0,
                  full::MlScore1);

DECLARE_SOA_TABLE(HfCandFlowInfos, "AOD", "HFCANDFLOWINFO",
                  full::M,
                  full::Pt,
                  full::MlScore0,
                  full::MlScore1,
                  full::ScalarProd,
                  full::Cent);
} // namespace o2::aod

enum DecayChannel { DplusToPiKPi = 0,
                    DsToKKPi,
                    DsToPiKK,
                    D0ToPiK,
                    D0ToKPi,
                    LcToPKPi,
                    LcToPiKP,
                    XicToPKPi,
                    XicToPiKP,
                    Xic0ToXiPi
};

struct HfTaskFlowCharmHadrons {
  Produces<o2::aod::HfCandMPtInfos> rowCandMassPtMl;
  Produces<o2::aod::HfCandFlowInfos> rowCandMassPtMlSpCent;

  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<int> qVecDetector{"qVecDetector", 3, "Detector for Q vector estimation (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3, TPC Pos: 4, TPC Neg: 5, TPC Tot: 6)"};
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3, FV0A: 4)"};
  Configurable<int> selectionFlag{"selectionFlag", 1, "Selection Flag for hadron (e.g. 1 for skimming, 3 for topo. and kine., 7 for PID)"};
  Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality accepted in SP/EP computation (not applied in resolution process)"};
  Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality accepted in SP/EP computation (not applied in resolution process)"};
  Configurable<bool> storeEP{"storeEP", false, "Flag to store EP-related axis"};
  Configurable<bool> storeMl{"storeMl", false, "Flag to store ML scores"};
  Configurable<bool> fillMassPtMlTree{"fillMassPtMlTree", false, "Flag to fill mass, pt and ML scores tree"};
  Configurable<bool> fillMassPtMlSpCentTree{"fillMassPtMlSpCentTree", false, "Flag to fill mass, pt, ML scores, SP and centrality tree"};
  Configurable<bool> fillSparse{"fillSparse", true, "Flag to fill sparse"};
  Configurable<float> downSampleFactor{"downSampleFactor", 1., "Fraction of candidates to keep in TTree"};
  Configurable<float> ptDownSampleMax{"ptDownSampleMax", 10., "Maximum pt for the application of the downsampling factor"};
  Configurable<bool> storeResoOccu{"storeResoOccu", false, "Flag to store Occupancy in resolution ThnSparse"};
  Configurable<bool> storeEpCosSin{"storeEpCosSin", false, "Flag to store cos and sin of EP angle in ThnSparse"};
  Configurable<bool> storeCandEta{"storeCandEta", false, "Flag to store candidates eta"};
  Configurable<bool> storeCandSign{"storeCandSign", false, "Flag to store candidates sign"};
  Configurable<int> occEstimator{"occEstimator", 0, "Occupancy estimation (0: None, 1: ITS, 2: FT0C)"};
  Configurable<bool> saveEpResoHisto{"saveEpResoHisto", false, "Flag to save event plane resolution histogram"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::vector<int>> classMl{"classMl", {0, 2}, "Indexes of BDT scores to be stored. Two indexes max."};

  EventPlaneHelper epHelper;
  HfEventSelection hfEvSel; // event selection and monitoring
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  SliceCache cache;

  using CandDsDataWMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfMlDsToKKPi>>;
  using CandDsData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi>>;
  using CandDplusDataWMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfMlDplusToPiKPi>>;
  using CandDplusData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi>>;
  using CandLcData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>;
  using CandLcDataWMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfMlLcToPKPi>>;
  using CandXicData = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi>>;
  using CandXicDataWMl = soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi, aod::HfMlXicToPKPi>>;
  using CandXic0Data = soa::Filtered<soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf>>;
  using CandXic0DataWMl = soa::Filtered<soa::Join<aod::HfCandToXiPiKf, aod::HfSelToXiPiKf, aod::HfMlToXiPi>>;
  using CandD0DataWMl = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfMlD0>>;
  using CandD0Data = soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>>;
  using CollsWithQvecs = soa::Join<aod::Collisions, aod::EvSels, aod::QvectorFT0Cs, aod::QvectorFT0As, aod::QvectorFT0Ms, aod::QvectorFV0As, aod::QvectorBPoss, aod::QvectorBNegs, aod::QvectorBTots, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>;
  using TracksWithExtra = soa::Join<aod::Tracks, aod::TracksExtra>;

  Filter filterSelectDsCandidates = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag || aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Filter filterSelectDplusCandidates = aod::hf_sel_candidate_dplus::isSelDplusToPiKPi >= selectionFlag;
  Filter filterSelectD0Candidates = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Filter filterSelectLcCandidates = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlag || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlag;
  Filter filterSelectXicCandidates = aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlag || aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlag;
  Filter filterSelectXic0Candidates = aod::hf_sel_toxipi::resultSelections == true;

  Partition<CandDsData> selectedDsToKKPi = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag;
  Partition<CandDsData> selectedDsToPiKK = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Partition<CandDsDataWMl> selectedDsToKKPiWMl = aod::hf_sel_candidate_ds::isSelDsToKKPi >= selectionFlag;
  Partition<CandDsDataWMl> selectedDsToPiKKWMl = aod::hf_sel_candidate_ds::isSelDsToPiKK >= selectionFlag;
  Partition<CandD0Data> selectedD0ToPiK = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
  Partition<CandD0Data> selectedD0ToKPi = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Partition<CandD0DataWMl> selectedD0ToPiKWMl = aod::hf_sel_candidate_d0::isSelD0 >= selectionFlag;
  Partition<CandD0DataWMl> selectedD0ToKPiWMl = aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlag;
  Partition<CandLcData> selectedLcToPKPi = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlag;
  Partition<CandLcData> selectedLcToPiKP = aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlag;
  Partition<CandLcDataWMl> selectedLcToPKPiWMl = aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlag;
  Partition<CandLcDataWMl> selectedLcToPiKPWMl = aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlag;
  Partition<CandXicData> selectedXicToPKPi = aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlag;
  Partition<CandXicData> selectedXicToPiKP = aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlag;
  Partition<CandXicDataWMl> selectedXicToPKPiWMl = aod::hf_sel_candidate_xic::isSelXicToPKPi >= selectionFlag;
  Partition<CandXicDataWMl> selectedXicToPiKPWMl = aod::hf_sel_candidate_xic::isSelXicToPiKP >= selectionFlag;
  Partition<CandXic0Data> selectedXic0 = aod::hf_sel_toxipi::resultSelections == true;
  Partition<CandXic0DataWMl> selectedXic0WMl = aod::hf_sel_toxipi::resultSelections == true;

  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {100, 1.78, 2.05}, ""};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {10, 0., 10.}, ""};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {10000, 0., 100.}, ""};
  ConfigurableAxis thnConfigAxisCosNPhi{"thnConfigAxisCosNPhi", {100, -1., 1.}, ""};
  ConfigurableAxis thnConfigAxisPsi{"thnConfigAxisPsi", {6000, 0, constants::math::TwoPI}, ""};
  ConfigurableAxis thnConfigAxisCosDeltaPhi{"thnConfigAxisCosDeltaPhi", {100, -1., 1.}, ""};
  ConfigurableAxis thnConfigAxisScalarProd{"thnConfigAxisScalarProd", {100, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisMlOne{"thnConfigAxisMlOne", {1000, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisMlTwo{"thnConfigAxisMlTwo", {1000, 0., 1.}, ""};
  ConfigurableAxis thnConfigAxisOccupancyITS{"thnConfigAxisOccupancyITS", {14, 0, 14000}, ""};
  ConfigurableAxis thnConfigAxisOccupancyFT0C{"thnConfigAxisOccupancyFT0C", {14, 0, 140000}, ""};
  ConfigurableAxis thnConfigAxisNoSameBunchPileup{"thnConfigAxisNoSameBunchPileup", {2, 0, 2}, ""};
  ConfigurableAxis thnConfigAxisOccupancy{"thnConfigAxisOccupancy", {2, 0, 2}, ""};
  ConfigurableAxis thnConfigAxisNoCollInTimeRangeNarrow{"thnConfigAxisNoCollInTimeRangeNarrow", {2, 0, 2}, ""};
  ConfigurableAxis thnConfigAxisNoCollInTimeRangeStandard{"thnConfigAxisNoCollInTimeRangeStandard", {2, 0, 2}, ""};
  ConfigurableAxis thnConfigAxisNoCollInRofStandard{"thnConfigAxisNoCollInRofStandard", {2, 0, 2}, ""};
  ConfigurableAxis thnConfigAxisResoFT0cFV0a{"thnConfigAxisResoFT0cFV0a", {160, -8, 8}, ""};
  ConfigurableAxis thnConfigAxisResoFT0cTPCtot{"thnConfigAxisResoFT0cTPCtot", {160, -8, 8}, ""};
  ConfigurableAxis thnConfigAxisResoFV0aTPCtot{"thnConfigAxisResoFV0aTPCtot", {160, -8, 8}, ""};
  ConfigurableAxis thnConfigAxisCandidateEta{"thnConfigAxisCandidateEta", {100, -5, 5}, ""};
  ConfigurableAxis thnConfigAxisSign{"thnConfigAxisSign", {6, -3.0, 3.0}, ""};

  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    if (storeResoOccu && occEstimator == 0) {
      LOGP(fatal, "Occupancy estimation must be enabled to store resolution THnSparse! Please check your configuration!");
    }
    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "Inv. mass (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality"};
    const AxisSpec thnAxisCosNPhi{thnConfigAxisCosNPhi, Form("cos(%d#varphi)", harmonic.value)};
    const AxisSpec thnAxisSinNPhi{thnConfigAxisCosNPhi, Form("sin(%d#varphi)", harmonic.value)};
    const AxisSpec thnAxisPsi{thnConfigAxisPsi, Form("#Psi_{%d}", harmonic.value)};
    const AxisSpec thnAxisCosDeltaPhi{thnConfigAxisCosDeltaPhi, Form("cos(%d(#varphi - #Psi_{sub}))", harmonic.value)};
    const AxisSpec thnAxisScalarProd{thnConfigAxisScalarProd, "SP"};
    const AxisSpec thnAxisMlOne{thnConfigAxisMlOne, "Bkg score"};
    const AxisSpec thnAxisMlTwo{thnConfigAxisMlTwo, "FD score"};
    const AxisSpec thnAxisOccupancyITS{thnConfigAxisOccupancyITS, "OccupancyITS"};
    const AxisSpec thnAxisOccupancyFT0C{thnConfigAxisOccupancyFT0C, "OccupancyFT0C"};
    const AxisSpec thnAxisCandEta{thnConfigAxisCandidateEta, "#eta"};
    const AxisSpec thnAxisSign{thnConfigAxisSign, "Sign"};
    const AxisSpec thnAxisNoSameBunchPileup{thnConfigAxisNoSameBunchPileup, "NoSameBunchPileup"};
    const AxisSpec thnAxisOccupancy{thnConfigAxisOccupancy, "Occupancy"};
    const AxisSpec thnAxisNoCollInTimeRangeNarrow{thnConfigAxisNoCollInTimeRangeNarrow, "NoCollInTimeRangeNarrow"};
    const AxisSpec thnAxisNoCollInTimeRangeStandard{thnConfigAxisNoCollInTimeRangeStandard, "NoCollInTimeRangeStandard"};
    const AxisSpec thnAxisNoCollInRofStandard{thnConfigAxisNoCollInRofStandard, "NoCollInRofStandard"};
    // TODO: currently only the Q vector of FT0c FV0a and TPCtot are considered
    const AxisSpec thnAxisResoFT0cFV0a{thnConfigAxisResoFT0cFV0a, "Q_{FT0c} #bullet Q_{FV0a}"};
    const AxisSpec thnAxisResoFT0cTPCtot{thnConfigAxisResoFT0cTPCtot, "Q_{FT0c} #bullet Q_{TPCtot}"};
    const AxisSpec thnAxisResoFV0aTPCtot{thnConfigAxisResoFV0aTPCtot, "Q_{FV0a} #bullet Q_{TPCtot}"};

    std::vector<AxisSpec> axes = {thnAxisInvMass, thnAxisPt, thnAxisCent, thnAxisScalarProd};
    if (storeEP) {
      axes.insert(axes.end(), {thnAxisCosNPhi, thnAxisSinNPhi, thnAxisCosDeltaPhi});
    }
    if (storeMl) {
      axes.insert(axes.end(), {thnAxisMlOne, thnAxisMlTwo});
    }
    if (storeCandEta) {
      axes.insert(axes.end(), {thnAxisCandEta});
    }
    if (storeCandSign) {
      axes.insert(axes.end(), {thnAxisSign});
    }
    if (occEstimator != 0) {
      if (occEstimator == 1) {
        axes.insert(axes.end(), {thnAxisOccupancyITS, thnAxisNoSameBunchPileup, thnAxisOccupancy,
                                 thnAxisNoCollInTimeRangeNarrow, thnAxisNoCollInTimeRangeStandard, thnAxisNoCollInRofStandard});
      } else {
        axes.insert(axes.end(), {thnAxisOccupancyFT0C, thnAxisNoSameBunchPileup, thnAxisOccupancy,
                                 thnAxisNoCollInTimeRangeNarrow, thnAxisNoCollInTimeRangeStandard, thnAxisNoCollInRofStandard});
      }
    }
    registry.add("hSparseFlowCharm", "THn for SP", HistType::kTHnSparseF, axes);

    if (occEstimator != 0) {
      registry.add("trackOccVsFT0COcc", "trackOccVsFT0COcc; trackOcc; FT0COcc", {HistType::kTH2F, {thnAxisOccupancyITS, thnAxisOccupancyFT0C}});
    }

    if (storeEpCosSin) {
      registry.add("ep/hSparseEp", "THn for Event Plane distirbution", {HistType::kTHnSparseF, {thnAxisCent, thnAxisPsi, thnAxisCosNPhi, thnAxisSinNPhi}});
    }

    if (doprocessResolution) { // enable resolution histograms only for resolution process
      registry.add("spReso/hSpResoFT0cFT0a", "hSpResoFT0cFT0a; centrality; Q_{FT0c} #bullet Q_{FT0a}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFT0cFV0a", "hSpResoFT0cFV0a; centrality; Q_{FT0c} #bullet Q_{FV0a}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFT0cTPCpos", "hSpResoFT0cTPCpos; centrality; Q_{FT0c} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFT0cTPCneg", "hSpResoFT0cTPCneg; centrality; Q_{FT0c} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFT0cTPCtot", "hSpResoFT0cTPCtot; centrality; Q_{FT0c} #bullet Q_{TPCtot}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFT0aFV0a", "hSpResoFT0aFV0a; centrality; Q_{FT0a} #bullet Q_{FV0a}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFT0aTPCpos", "hSpResoFT0aTPCpos; centrality; Q_{FT0a} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFT0aTPCneg", "hSpResoFT0aTPCneg; centrality; Q_{FT0a} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFT0aTPCtot", "hSpResoFT0aTPCtot; centrality; Q_{FT0m} #bullet Q_{TPCtot}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFT0mFV0a", "hSpResoFT0mFV0a; centrality; Q_{FT0m} #bullet Q_{FV0a}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFT0mTPCpos", "hSpResoFT0mTPCpos; centrality; Q_{FT0m} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFT0mTPCneg", "hSpResoFT0mTPCneg; centrality; Q_{FT0m} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFT0mTPCtot", "hSpResoFT0mTPCtot; centrality; Q_{FV0a} #bullet Q_{TPCtot}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFV0aTPCpos", "hSpResoFV0aTPCpos; centrality; Q_{FV0a} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFV0aTPCneg", "hSpResoFV0aTPCneg; centrality; Q_{FV0a} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoFV0aTPCtot", "hSpResoFV0aTPCtot; centrality; Q_{FV0a} #bullet Q_{TPCtot}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
      registry.add("spReso/hSpResoTPCposTPCneg", "hSpResoTPCposTPCneg; centrality; Q_{TPCpos} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});

      if (saveEpResoHisto) {
        registry.add("epReso/hEpResoFT0cFT0a", "hEpResoFT0cFT0a; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFT0cFV0a", "hEpResoFT0cFV0a; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFT0cTPCpos", "hEpResoFT0cTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFT0cTPCneg", "hEpResoFT0cTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFT0cTPCtot", "hEpResoFT0cTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFT0aFV0a", "hEpResoFT0aFV0a; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFT0aTPCpos", "hEpResoFT0aTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFT0aTPCneg", "hEpResoFT0aTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFT0aTPCtot", "hEpResoFT0aTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFT0mFV0a", "hEpResoFT0mFV0a; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFT0mTPCpos", "hEpResoFT0mTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFT0mTPCneg", "hEpResoFT0mTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFT0mTPCtot", "hEpResoFT0mTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFV0aTPCpos", "hEpResoFV0aTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFV0aTPCneg", "hEpResoFV0aTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoFV0aTPCtot", "hEpResoFV0aTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
        registry.add("epReso/hEpResoTPCposTPCneg", "hEpResoTPCposTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTH2F, {thnAxisCent, thnAxisCosNPhi}});
      }

      if (storeResoOccu) {
        std::vector<AxisSpec> axesReso = {thnAxisCent, thnAxisResoFT0cFV0a, thnAxisResoFT0cTPCtot, thnAxisResoFV0aTPCtot};
        if (occEstimator == 1) {
          axesReso.insert(axesReso.end(), {thnAxisOccupancyITS, thnAxisNoSameBunchPileup, thnAxisOccupancy,
                                           thnAxisNoCollInTimeRangeNarrow, thnAxisNoCollInTimeRangeStandard, thnAxisNoCollInRofStandard});
        } else {
          axesReso.insert(axesReso.end(), {thnAxisOccupancyFT0C, thnAxisNoSameBunchPileup, thnAxisOccupancy,
                                           thnAxisNoCollInTimeRangeNarrow, thnAxisNoCollInTimeRangeStandard, thnAxisNoCollInRofStandard});
        }
        registry.add("spReso/hSparseReso", "THn for resolution with occupancy", HistType::kTHnSparseF, axesReso);
      }

      hfEvSel.init(registry);
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
    }
  }; // end init

  /// Compute the Q vector for the candidate's tracks
  /// \param cand is the candidate
  /// \param tracksQx is the X component of the Q vector for the tracks
  /// \param tracksQy is the Y component of the Q vector for the tracks
  /// \param channel is the decay channel
  template <DecayChannel Channel, typename T1>
  void getQvecDtracks(T1 const& cand,
                      std::vector<float>& tracksQx,
                      std::vector<float>& tracksQy,
                      const float amplQVec,
                      QvecEstimator qVecDetector)
  {
    auto addProngIfInSubevent = [&](float px, float py, float pz) {
      const std::array<float, 3> pVec{px, py, pz};
      const float eta = RecoDecay::eta(pVec);

      // only subtract daughters that actually contributed to THIS subevent Q
      if (qVecDetector == QvecEstimator::TPCPos && eta <= 0.f) {
        return;
      }
      if (qVecDetector == QvecEstimator::TPCNeg && eta >= 0.f) {
        return;
      }
      // for TPCTot: no early return, all prongs contribute

      const float pt = RecoDecay::pt(pVec); // or std::hypot(px, py);
      const float phi = std::atan2(py, px);

      // xQVec,yQVec are normalized with amplQVec, so the daughters' contribution must use the SAME normalization.
      tracksQx.push_back(std::cos(harmonic * phi) * pt / amplQVec);
      tracksQy.push_back(std::sin(harmonic * phi) * pt / amplQVec);
    };

    addProngIfInSubevent(cand.pxProng0(), cand.pyProng0(), cand.pzProng0());
    addProngIfInSubevent(cand.pxProng1(), cand.pyProng1(), cand.pzProng1());

    // 3-prong channels
    if constexpr (Channel != DecayChannel::D0ToPiK && Channel != DecayChannel::D0ToKPi) {
      addProngIfInSubevent(cand.pxProng2(), cand.pyProng2(), cand.pzProng2());
    }
  }

  /// Compute the Q vector for the candidate's tracks
  /// \param cand is the candidate
  /// \param tracksQx is the X component of the Q vector for the tracks
  /// \param tracksQy is the Y component of the Q vector for the tracks
  /// \param channel is the decay channel
  template <typename T1>
  void getQvecXic0Tracks(const T1& cand,
                         std::vector<float>& tracksQx,
                         std::vector<float>& tracksQy,
                         float amplQVec,
                         QvecEstimator qVecDetector)
  {
    auto addProngIfInSubevent = [&](float px, float py, float pz) {
      const std::array<float, 3> pVec{px, py, pz};
      const float eta = RecoDecay::eta(pVec);

      if (qVecDetector == QvecEstimator::TPCPos && eta <= 0.f) {
        return;
      }
      if (qVecDetector == QvecEstimator::TPCNeg && eta >= 0.f) {
        return;
      }

      const float pt = RecoDecay::pt(pVec);
      const float phi = std::atan2(py, px);

      tracksQx.push_back(std::cos(harmonic * phi) * pt / amplQVec);
      tracksQy.push_back(std::sin(harmonic * phi) * pt / amplQVec);
    };

    addProngIfInSubevent(cand.pxPosV0Dau(), cand.pyPosV0Dau(), cand.pzPosV0Dau());
    addProngIfInSubevent(cand.pxNegV0Dau(), cand.pyNegV0Dau(), cand.pzNegV0Dau());
    addProngIfInSubevent(cand.pxBachFromCasc(), cand.pyBachFromCasc(), cand.pzBachFromCasc());
    addProngIfInSubevent(cand.pxBachFromCharmBaryon(), cand.pyBachFromCharmBaryon(), cand.pzBachFromCharmBaryon());
  }

  /// Get the event selection flags
  /// \param hfevselflag is the event selection flag
  std::vector<int> getEventSelectionFlags(const o2::hf_evsel::HfCollisionRejectionMask hfevselflag)
  {
    return {
      TESTBIT(hfevselflag, o2::hf_evsel::EventRejection::NoSameBunchPileup),
      TESTBIT(hfevselflag, o2::hf_evsel::EventRejection::Occupancy),
      TESTBIT(hfevselflag, o2::hf_evsel::EventRejection::NoCollInTimeRangeNarrow),
      TESTBIT(hfevselflag, o2::hf_evsel::EventRejection::NoCollInTimeRangeStandard),
      TESTBIT(hfevselflag, o2::hf_evsel::EventRejection::NoCollInRofStandard)};
  }

  /// Fill THnSparse
  /// \param mass is the invariant mass of the candidate
  /// \param pt is the transverse momentum of the candidate
  /// \param eta is the pseudorapidityof the candidate
  /// \param sign is the particle charge sign of the candidate
  /// \param cent is the centrality of the collision
  /// \param cosNPhi is the cosine of the n*phi angle
  /// \param sinNPhi is the sine of the n*phi angle
  /// \param cosDeltaPhi is the cosine of the n*(phi - evtPl) angle
  /// \param sp is the scalar product
  /// \param outputMl are the ML scores
  /// \param occupancy is the occupancy of the collision using the track estimator
  /// \param hfevselflag flag of the collision associated to utilsEvSelHf.h
  void fillThn(const float mass,
               const float pt,
               const float eta,
               const float sign,
               const float cent,
               const float cosNPhi,
               const float sinNPhi,
               const float cosDeltaPhi,
               const float sp,
               const std::vector<float>& outputMl,
               const float occupancy,
               const o2::hf_evsel::HfCollisionRejectionMask hfevselflag)
  {
    auto hSparse = registry.get<THnSparse>(HIST("hSparseFlowCharm"));
    const int ndim = hSparse->GetNdimensions();

    std::vector<double> values;
    values.reserve(ndim);

    values.push_back(mass);
    values.push_back(pt);
    values.push_back(cent);
    values.push_back(sp);

    if (storeEP) {
      values.push_back(cosNPhi);
      values.push_back(sinNPhi);
      values.push_back(cosDeltaPhi);
    }

    if (storeMl) {
      values.push_back(outputMl[0]);
      values.push_back(outputMl[1]);
    }

    if (storeCandEta) {
      values.push_back(eta);
    }
    if (storeCandSign) {
      values.push_back(sign);
    }
    if (occEstimator != 0) {
      auto evtSelFlags = getEventSelectionFlags(hfevselflag);
      values.push_back(occupancy);
      values.push_back(evtSelFlags[0]);
      values.push_back(evtSelFlags[1]);
      values.push_back(evtSelFlags[2]);
      values.push_back(evtSelFlags[3]);
      values.push_back(evtSelFlags[4]);
    }

    if (static_cast<int>(values.size()) != ndim) {
      LOGF(fatal,
           "hSparseFlowCharm: number of filled dimensions (%d) "
           "does not match THnSparse dimensionality (%d).",
           static_cast<int>(values.size()), ndim);
    }

    hSparse->Fill(values.data());
  }

  /// Check if the collision is selected
  /// \param collision is the collision with the Q vector information
  /// \param bc is the bunch crossing with timestamp information
  /// \param centrality is the collision centrality
  /// \return true if the collision is selected, false otherwise
  template <o2::hf_centrality::CentralityEstimator CentEstimator>
  bool isCollSelected(CollsWithQvecs::iterator const& collision,
                      aod::BCsWithTimestamps const&,
                      float& centrality)
  {
    const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
    const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
    centrality = o2::hf_centrality::getCentralityColl(collision, CentEstimator);

    /// monitor the satisfied event selections
    hfEvSel.fillHistograms(collision, rejectionMask, centrality, occupancy);
    registry.fill(HIST("trackOccVsFT0COcc"), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange());
    return rejectionMask == 0;
  }

  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param candidates are the selected candidates
  template <DecayChannel Channel, typename T1, typename Trk>
  void runFlowAnalysis(CollsWithQvecs::iterator const& collision,
                       T1 const& candidates,
                       Trk const& /*tracks*/)
  {
    float cent = o2::hf_centrality::getCentralityColl(collision, centEstimator);
    if (cent < centralityMin || cent > centralityMax) {
      return;
    }
    float occupancy = 0.;
    o2::hf_evsel::HfCollisionRejectionMask hfevflag{};
    if (occEstimator != 0) {
      occupancy = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
      registry.fill(HIST("trackOccVsFT0COcc"), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange());
      hfevflag = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, cent, ccdb, registry);
    }

    std::array<float, 3> qVecs = getQvec(collision, qVecDetector.value);
    float xQVec = qVecs[0];
    float yQVec = qVecs[1];
    float const amplQVec = qVecs[2];
    float const evtPl = epHelper.GetEventPlane(xQVec, yQVec, harmonic);
    for (const auto& candidate : candidates) {
      float massCand = 0.;
      float signCand = 0.;
      std::vector<float> outputMl = {-999., -999.};
      auto trackprong0 = candidate.template prong0_as<Trk>();
      if constexpr (std::is_same_v<T1, CandDsData> || std::is_same_v<T1, CandDsDataWMl>) {
        signCand = trackprong0.sign();
        switch (Channel) {
          case DecayChannel::DsToKKPi:
            massCand = HfHelper::invMassDsToKKPi(candidate);
            if constexpr (std::is_same_v<T1, CandDsDataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
                outputMl[iclass] = candidate.mlProbDsToKKPi()[classMl->at(iclass)];
              }
            }
            break;
          case DecayChannel::DsToPiKK:
            massCand = HfHelper::invMassDsToPiKK(candidate);
            if constexpr (std::is_same_v<T1, CandDsDataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
                outputMl[iclass] = candidate.mlProbDsToPiKK()[classMl->at(iclass)];
              }
            }
            break;
          default:
            break;
        }
      } else if constexpr (std::is_same_v<T1, CandDplusData> || std::is_same_v<T1, CandDplusDataWMl>) {
        massCand = HfHelper::invMassDplusToPiKPi(candidate);
        signCand = trackprong0.sign();
        if constexpr (std::is_same_v<T1, CandDplusDataWMl>) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
          }
        }
      } else if constexpr (std::is_same_v<T1, CandD0Data> || std::is_same_v<T1, CandD0DataWMl>) {
        switch (Channel) {
          case DecayChannel::D0ToPiK:
            signCand = candidate.isSelD0bar() ? 3 : 1; // 3: reflected D0bar, 1: pure D0 excluding reflected D0bar
            massCand = HfHelper::invMassD0ToPiK(candidate);
            if constexpr (std::is_same_v<T1, CandD0DataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
                outputMl[iclass] = candidate.mlProbD0()[classMl->at(iclass)];
              }
            }
            break;
          case DecayChannel::D0ToKPi:
            massCand = HfHelper::invMassD0barToKPi(candidate);
            signCand = candidate.isSelD0() ? 3 : 2; // 3: reflected D0, 2: pure D0bar excluding reflected D0
            if constexpr (std::is_same_v<T1, CandD0DataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
                outputMl[iclass] = candidate.mlProbD0bar()[classMl->at(iclass)];
              }
            }
            break;
          default:
            break;
        }
      } else if constexpr (std::is_same_v<T1, CandLcData> || std::is_same_v<T1, CandLcDataWMl>) {
        signCand = trackprong0.sign();
        switch (Channel) {
          case DecayChannel::LcToPKPi:
            massCand = HfHelper::invMassLcToPKPi(candidate);
            if constexpr (std::is_same_v<T1, CandLcDataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
                outputMl[iclass] = candidate.mlProbLcToPKPi()[classMl->at(iclass)];
              }
            }
            break;
          case DecayChannel::LcToPiKP:
            massCand = HfHelper::invMassLcToPiKP(candidate);
            if constexpr (std::is_same_v<T1, CandLcDataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
                outputMl[iclass] = candidate.mlProbLcToPiKP()[classMl->at(iclass)];
              }
            }
            break;
          default:
            break;
        }
      } else if constexpr (std::is_same_v<T1, CandXicData> || std::is_same_v<T1, CandXicDataWMl>) {
        signCand = trackprong0.sign();
        switch (Channel) {
          case DecayChannel::XicToPKPi:
            massCand = HfHelper::invMassXicToPKPi(candidate);
            if constexpr (std::is_same_v<T1, CandXicDataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
                outputMl[iclass] = candidate.mlProbXicToPKPi()[classMl->at(iclass)];
              }
            }
            break;
          case DecayChannel::XicToPiKP:
            massCand = HfHelper::invMassXicToPiKP(candidate);
            if constexpr (std::is_same_v<T1, CandXicDataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
                outputMl[iclass] = candidate.mlProbXicToPiKP()[classMl->at(iclass)];
              }
            }
            break;
          default:
            break;
        }
      } else if constexpr (std::is_same_v<T1, CandXic0Data> || std::is_same_v<T1, CandXic0DataWMl>) {
        massCand = candidate.invMassCharmBaryon();
        signCand = static_cast<float>(candidate.signDecay());
        if constexpr (std::is_same_v<T1, CandXic0DataWMl>) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbToXiPi()[classMl->at(iclass)];
          }
        }
      }

      if constexpr (std::is_same_v<T1, CandXic0Data> || std::is_same_v<T1, CandXic0DataWMl>) {
        signCand = candidate.signDecay();
      } else {
        auto trackprong0 = candidate.template prong0_as<Trk>();
        signCand = trackprong0.sign();
      }

      float ptCand = 0.;
      float phiCand = 0.;
      float etaCand = 0.;

      if constexpr (std::is_same_v<T1, CandXic0Data> || std::is_same_v<T1, CandXic0DataWMl>) {
        ptCand = RecoDecay::sqrtSumOfSquares(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
        phiCand = std::atan2(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
        etaCand = candidate.etaBachFromCharmBaryon();
      } else {
        ptCand = candidate.pt();
        phiCand = candidate.phi();
        etaCand = candidate.eta();
      }

      // If TPC is used for the SP estimation, the tracks of the hadron candidate must be removed from the corresponding TPC Q vector to avoid self-correlations
      if (qVecDetector == QvecEstimator::TPCNeg ||
          qVecDetector == QvecEstimator::TPCPos ||
          qVecDetector == QvecEstimator::TPCTot) {

        std::vector<float> tracksQx;
        std::vector<float> tracksQy;

        // IMPORTANT: use the ORIGINAL amplitude used to build this Q-vector
        if constexpr (std::is_same_v<T1, CandXic0Data> || std::is_same_v<T1, CandXic0DataWMl>) {
          getQvecXic0Tracks(candidate, tracksQx, tracksQy, amplQVec, static_cast<QvecEstimator>(qVecDetector.value));
        } else {
          getQvecDtracks<Channel>(candidate, tracksQx, tracksQy, amplQVec, static_cast<QvecEstimator>(qVecDetector.value));
        }

        // subtract daughters' contribution from the (normalized) Q-vector
        for (std::size_t iTrack = 0; iTrack < tracksQx.size(); ++iTrack) {
          xQVec -= tracksQx[iTrack];
          yQVec -= tracksQy[iTrack];
        }
      }

      float const cosNPhi = std::cos(harmonic * phiCand);
      float const sinNPhi = std::sin(harmonic * phiCand);
      float const scalprodCand = cosNPhi * xQVec + sinNPhi * yQVec;
      float const cosDeltaPhi = std::cos(harmonic * (phiCand - evtPl));

      if (fillMassPtMlTree || fillMassPtMlSpCentTree) {
        if (downSampleFactor < 1.) {
          float const pseudoRndm = ptCand * 1000. - static_cast<int64_t>(ptCand * 1000);
          if (ptCand < ptDownSampleMax && pseudoRndm >= downSampleFactor) {
            continue;
          }
        }
        if (fillMassPtMlTree) {
          rowCandMassPtMl(massCand, ptCand, outputMl[0], outputMl[1]);
        }
        if (fillMassPtMlSpCentTree) {
          rowCandMassPtMlSpCent(massCand, ptCand, outputMl[0], outputMl[1], scalprodCand, cent);
        }
      }
      if (fillSparse) {
        fillThn(massCand, ptCand, etaCand, signCand, cent, cosNPhi, sinNPhi, cosDeltaPhi, scalprodCand, outputMl, occupancy, hfevflag);
      }
    }
  }

  // Ds with ML
  void processDsMl(CollsWithQvecs::iterator const& collision,
                   CandDsDataWMl const&,
                   TracksWithExtra const& tracks)
  {
    auto candsDsToKKPiWMl = selectedDsToKKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsDsToPiKKWMl = selectedDsToPiKKWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::DsToKKPi>(collision, candsDsToKKPiWMl, tracks);
    runFlowAnalysis<DecayChannel::DsToPiKK>(collision, candsDsToPiKKWMl, tracks);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDsMl, "Process Ds candidates with ML", false);

  // Ds with rectangular cuts
  void processDs(CollsWithQvecs::iterator const& collision,
                 CandDsData const& /*candidatesDs*/,
                 TracksWithExtra const& tracks)
  {
    auto candsDsToKKPi = selectedDsToKKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsDsToPiKK = selectedDsToPiKK->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::DsToKKPi>(collision, candsDsToKKPi, tracks);
    runFlowAnalysis<DecayChannel::DsToPiKK>(collision, candsDsToPiKK, tracks);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDs, "Process Ds candidates", false);

  // Dplus with ML
  void processDplusMl(CollsWithQvecs::iterator const& collision,
                      CandDplusDataWMl const& candidatesDplus,
                      TracksWithExtra const& tracks)
  {
    runFlowAnalysis<DecayChannel::DplusToPiKPi>(collision, candidatesDplus, tracks);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDplusMl, "Process Dplus candidates with ML", false);

  // Dplus with rectangular cuts
  void processDplus(CollsWithQvecs::iterator const& collision,
                    CandDplusData const& candidatesDplus,
                    TracksWithExtra const& tracks)
  {
    runFlowAnalysis<DecayChannel::DplusToPiKPi>(collision, candidatesDplus, tracks);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDplus, "Process Dplus candidates", true);

  // D0 with ML
  void processD0Ml(CollsWithQvecs::iterator const& collision,
                   CandD0DataWMl const& /*candidatesD0*/,
                   TracksWithExtra const& tracks)
  {
    auto candsD0ToPiKWMl = selectedD0ToPiKWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsD0ToKPiWMl = selectedD0ToKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::D0ToPiK>(collision, candsD0ToPiKWMl, tracks);
    runFlowAnalysis<DecayChannel::D0ToKPi>(collision, candsD0ToKPiWMl, tracks);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processD0Ml, "Process D0 candidates with ML", false);

  // D0 with rectangular cuts
  void processD0(CollsWithQvecs::iterator const& collision,
                 CandD0Data const& /*candidatesD0*/,
                 TracksWithExtra const& tracks)
  {
    auto candsD0ToPiK = selectedD0ToPiK->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsD0ToKPi = selectedD0ToKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::D0ToPiK>(collision, candsD0ToPiK, tracks);
    runFlowAnalysis<DecayChannel::D0ToKPi>(collision, candsD0ToKPi, tracks);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processD0, "Process D0 candidates", false);

  // Lc with ML
  void processLcMl(CollsWithQvecs::iterator const& collision,
                   CandLcDataWMl const& /*candidatesLc*/,
                   TracksWithExtra const& tracks)
  {
    auto candsLcToPKPiWMl = selectedLcToPKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsLcToPiKPWMl = selectedLcToPiKPWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::LcToPKPi>(collision, candsLcToPKPiWMl, tracks);
    runFlowAnalysis<DecayChannel::LcToPiKP>(collision, candsLcToPiKPWMl, tracks);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processLcMl, "Process Lc candidates with ML", false);

  // Lc with rectangular cuts
  void processLc(CollsWithQvecs::iterator const& collision,
                 CandLcData const& /*candidatesLc*/,
                 TracksWithExtra const& tracks)
  {
    auto candsLcToPKPi = selectedLcToPKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsLcToPiKP = selectedLcToPiKP->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::LcToPKPi>(collision, candsLcToPKPi, tracks);
    runFlowAnalysis<DecayChannel::LcToPiKP>(collision, candsLcToPiKP, tracks);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processLc, "Process Lc candidates", false);

  // Xic with ML
  void processXicMl(CollsWithQvecs::iterator const& collision,
                    CandXicDataWMl const& /*candidatesXic*/,
                    TracksWithExtra const& tracks)
  {
    auto candsXicToPKPiWMl = selectedXicToPKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsXicToPiKPWMl = selectedXicToPiKPWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::XicToPKPi>(collision, candsXicToPKPiWMl, tracks);
    runFlowAnalysis<DecayChannel::XicToPiKP>(collision, candsXicToPiKPWMl, tracks);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processXicMl, "Process Xic candidates with ML", false);

  // Xic with rectangular cuts
  void processXic(CollsWithQvecs::iterator const& collision,
                  CandXicData const& /*candidatesXic*/,
                  TracksWithExtra const& tracks)
  {
    auto candsXicToPKPi = selectedXicToPKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsXicToPiKP = selectedXicToPiKP->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::XicToPKPi>(collision, candsXicToPKPi, tracks);
    runFlowAnalysis<DecayChannel::XicToPiKP>(collision, candsXicToPiKP, tracks);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processXic, "Process Xic candidates", false);

  // Xic0 with ML
  void processXic0Ml(CollsWithQvecs::iterator const& collision,
                     CandXic0DataWMl const& /*candidatesXic0*/,
                     TracksWithExtra const& tracks)
  {
    auto candsXic0WMl = selectedXic0WMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::Xic0ToXiPi>(collision, candsXic0WMl, tracks);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processXic0Ml, "Process Xic0 candidates with ML", false);

  // Xic0
  void processXic0(CollsWithQvecs::iterator const& collision,
                   CandXic0Data const& /*candidatesXic0*/,
                   TracksWithExtra const& tracks)
  {
    auto candsXic0 = selectedXic0->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::Xic0ToXiPi>(collision, candsXic0, tracks);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processXic0, "Process Xic0 candidates", false);

  // Resolution
  void processResolution(CollsWithQvecs::iterator const& collision,
                         aod::BCsWithTimestamps const& bcs)
  {
    float centrality{-1.f};
    float const xQVecFT0a = collision.qvecFT0ARe();
    float const yQVecFT0a = collision.qvecFT0AIm();
    float const xQVecFT0c = collision.qvecFT0CRe();
    float const yQVecFT0c = collision.qvecFT0CIm();
    float const xQVecFT0m = collision.qvecFT0MRe();
    float const yQVecFT0m = collision.qvecFT0MIm();
    float const xQVecFV0a = collision.qvecFV0ARe();
    float const yQVecFV0a = collision.qvecFV0AIm();
    float const xQVecBPos = collision.qvecBPosRe();
    float const yQVecBPos = collision.qvecBPosIm();
    float const xQVecBNeg = collision.qvecBNegRe();
    float const yQVecBNeg = collision.qvecBNegIm();
    float const xQVecBTot = collision.qvecBTotRe();
    float const yQVecBTot = collision.qvecBTotIm();

    centrality = o2::hf_centrality::getCentralityColl(collision, o2::hf_centrality::CentralityEstimator::FT0C);
    if (storeResoOccu) {
      const auto occupancy = o2::hf_occupancy::getOccupancyColl(collision, occEstimator);
      registry.fill(HIST("trackOccVsFT0COcc"), collision.trackOccupancyInTimeRange(), collision.ft0cOccupancyInTimeRange());
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, o2::hf_centrality::CentralityEstimator::None, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      std::vector<int> evtSelFlags = getEventSelectionFlags(rejectionMask);
      registry.fill(HIST("spReso/hSparseReso"), centrality, xQVecFT0c * xQVecFV0a + yQVecFT0c * yQVecFV0a,
                    xQVecFT0c * xQVecBTot + yQVecFT0c * yQVecBTot,
                    xQVecFV0a * xQVecBTot + yQVecFV0a * yQVecBTot,
                    occupancy, evtSelFlags[0], evtSelFlags[1], evtSelFlags[2], evtSelFlags[3], evtSelFlags[4]);
    }

    if (!isCollSelected<o2::hf_centrality::CentralityEstimator::FT0C>(collision, bcs, centrality)) {
      // no selection on the centrality is applied, but on event selection flags
      return;
    }

    registry.fill(HIST("spReso/hSpResoFT0cFT0a"), centrality, xQVecFT0c * xQVecFT0a + yQVecFT0c * yQVecFT0a);
    registry.fill(HIST("spReso/hSpResoFT0cFV0a"), centrality, xQVecFT0c * xQVecFV0a + yQVecFT0c * yQVecFV0a);
    registry.fill(HIST("spReso/hSpResoFT0cTPCpos"), centrality, xQVecFT0c * xQVecBPos + yQVecFT0c * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFT0cTPCneg"), centrality, xQVecFT0c * xQVecBNeg + yQVecFT0c * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFT0cTPCtot"), centrality, xQVecFT0c * xQVecBTot + yQVecFT0c * yQVecBTot);
    registry.fill(HIST("spReso/hSpResoFT0aFV0a"), centrality, xQVecFT0a * xQVecFV0a + yQVecFT0a * yQVecFV0a);
    registry.fill(HIST("spReso/hSpResoFT0aTPCpos"), centrality, xQVecFT0a * xQVecBPos + yQVecFT0a * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFT0aTPCneg"), centrality, xQVecFT0a * xQVecBNeg + yQVecFT0a * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFT0aTPCtot"), centrality, xQVecFT0a * xQVecBTot + yQVecFT0a * yQVecBTot);
    registry.fill(HIST("spReso/hSpResoFT0mFV0a"), centrality, xQVecFT0m * xQVecFV0a + yQVecFT0m * yQVecFV0a);
    registry.fill(HIST("spReso/hSpResoFT0mTPCpos"), centrality, xQVecFT0m * xQVecBPos + yQVecFT0m * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFT0mTPCneg"), centrality, xQVecFT0m * xQVecBNeg + yQVecFT0m * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFT0mTPCtot"), centrality, xQVecFT0m * xQVecBTot + yQVecFT0m * yQVecBTot);
    registry.fill(HIST("spReso/hSpResoFV0aTPCpos"), centrality, xQVecFV0a * xQVecBPos + yQVecFV0a * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFV0aTPCneg"), centrality, xQVecFV0a * xQVecBNeg + yQVecFV0a * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFV0aTPCtot"), centrality, xQVecFV0a * xQVecBTot + yQVecFV0a * yQVecBTot);
    registry.fill(HIST("spReso/hSpResoTPCposTPCneg"), centrality, xQVecBPos * xQVecBNeg + yQVecBPos * yQVecBNeg);

    if (saveEpResoHisto) {
      float const epFT0a = epHelper.GetEventPlane(xQVecFT0a, yQVecFT0a, harmonic);
      float const epFT0c = epHelper.GetEventPlane(xQVecFT0c, yQVecFT0c, harmonic);
      float const epFT0m = epHelper.GetEventPlane(xQVecFT0m, yQVecFT0m, harmonic);
      float const epFV0a = epHelper.GetEventPlane(xQVecFV0a, yQVecFV0a, harmonic);
      float const epBPoss = epHelper.GetEventPlane(xQVecBPos, yQVecBPos, harmonic);
      float const epBNegs = epHelper.GetEventPlane(xQVecBNeg, yQVecBNeg, harmonic);
      float const epBTots = epHelper.GetEventPlane(xQVecBTot, yQVecBTot, harmonic);

      registry.fill(HIST("epReso/hEpResoFT0cFT0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epFT0a, harmonic)));
      registry.fill(HIST("epReso/hEpResoFT0cFV0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epFV0a, harmonic)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBPoss, harmonic)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBNegs, harmonic)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBTots, harmonic)));
      registry.fill(HIST("epReso/hEpResoFT0aFV0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epFV0a, harmonic)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBPoss, harmonic)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBNegs, harmonic)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBTots, harmonic)));
      registry.fill(HIST("epReso/hEpResoFT0mFV0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epFV0a, harmonic)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBPoss, harmonic)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBNegs, harmonic)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBTots, harmonic)));
      registry.fill(HIST("epReso/hEpResoFV0aTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epBPoss, harmonic)));
      registry.fill(HIST("epReso/hEpResoFV0aTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epBNegs, harmonic)));
      registry.fill(HIST("epReso/hEpResoFV0aTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epBTots, harmonic)));
      registry.fill(HIST("epReso/hEpResoTPCposTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epBPoss, epBNegs, harmonic)));
    }

    if (storeEpCosSin) {
      registry.fill(HIST("ep/hSparseEp"), centrality,
                    epHelper.GetEventPlane(xQVecFT0c, yQVecFT0c, harmonic),
                    std::cos(harmonic * epHelper.GetEventPlane(xQVecFT0c, yQVecFT0c, harmonic)),
                    std::sin(harmonic * epHelper.GetEventPlane(xQVecFT0c, yQVecFT0c, harmonic)));
    }
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processResolution, "Process resolution", false);

}; // End struct HfTaskFlowCharmHadrons

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskFlowCharmHadrons>(cfgc)};
}
