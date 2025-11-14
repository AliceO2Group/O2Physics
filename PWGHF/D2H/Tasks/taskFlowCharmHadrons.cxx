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

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(M, m, float);   //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float); //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(MlScore0, mlScore0, float); //! ML score of the first configured index
DECLARE_SOA_COLUMN(MlScore1, mlScore1, float); //! ML score of the second configured index
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

enum QvecEstimator { FV0A = 0,
                     FT0M,
                     FT0A,
                     FT0C,
                     TPCPos,
                     TPCNeg,
                     TPCTot };

struct HfTaskFlowCharmHadrons {
  Produces<o2::aod::HfCandMPtInfos> rowCandMassPtMl;
  Produces<o2::aod::HfCandFlowInfos> rowCandMassPtMlSpCent;

  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<int> qvecDetector{"qvecDetector", 3, "Detector for Q vector estimation (FV0A: 0, FT0M: 1, FT0A: 2, FT0C: 3, TPC Pos: 4, TPC Neg: 5, TPC Tot: 6)"};
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

      hfEvSel.addHistograms(registry); // collision monitoring
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
                      const float ampl)
  {
    // TODO: add possibility to consider different weights for the tracks, at the moment only pT is considered;
    float const pXTrack0 = cand.pxProng0();
    float const pYTrack0 = cand.pyProng0();
    float const pTTrack0 = cand.ptProng0();
    float const phiTrack0 = std::atan2(pYTrack0, pXTrack0);
    float const pXTrack1 = cand.pxProng1();
    float const pYTrack1 = cand.pyProng1();
    float const pTTrack1 = cand.ptProng1();
    float const phiTrack1 = std::atan2(pYTrack1, pXTrack1);

    tracksQx.push_back(std::cos(harmonic * phiTrack0) * pTTrack0 / ampl);
    tracksQy.push_back(std::sin(harmonic * phiTrack0) * pTTrack0 / ampl);
    tracksQx.push_back(std::cos(harmonic * phiTrack1) * pTTrack1 / ampl);
    tracksQy.push_back(std::sin(harmonic * phiTrack1) * pTTrack1 / ampl);

    if constexpr (Channel != DecayChannel::D0ToPiK && Channel != DecayChannel::D0ToKPi) {
      float const pXTrack2 = cand.pxProng2();
      float const pYTrack2 = cand.pyProng2();
      float const pTTrack2 = cand.ptProng2();
      float const phiTrack2 = std::atan2(pYTrack2, pXTrack2);
      tracksQx.push_back(std::cos(harmonic * phiTrack2) * pTTrack2 / ampl);
      tracksQy.push_back(std::sin(harmonic * phiTrack2) * pTTrack2 / ampl);
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
                         float ampl)
  {
    // add possibility to consider different weights for the tracks, at the moment only pT is considered;
    float const pXTrack0 = cand.pxPosV0Dau();
    float const pYTrack0 = cand.pyPosV0Dau();
    float const pTTrack0 = std::hypot(pXTrack0, pYTrack0);
    float const phiTrack0 = std::atan2(pXTrack0, pYTrack0);
    float const pXTrack1 = cand.pxNegV0Dau();
    float const pYTrack1 = cand.pyNegV0Dau();
    float const pTTrack1 = std::hypot(pXTrack1, pYTrack1);
    float const phiTrack1 = std::atan2(pXTrack1, pYTrack1);
    float const pYTrack2 = cand.pxBachFromCasc();
    float const pXTrack2 = cand.pyBachFromCasc();
    float const pTTrack2 = std::hypot(pXTrack2, pYTrack2);
    float const phiTrack2 = std::atan2(pXTrack2, pYTrack2);
    float const pXTrack3 = cand.pxBachFromCharmBaryon();
    float const pYTrack3 = cand.pyBachFromCharmBaryon();
    float const pTTrack3 = std::hypot(pXTrack3, pYTrack3);
    float const phiTrack3 = std::atan2(pXTrack3, pYTrack3);

    tracksQx.push_back(std::cos(harmonic * phiTrack0) * pTTrack0 / ampl);
    tracksQy.push_back(std::sin(harmonic * phiTrack0) * pTTrack0 / ampl);
    tracksQx.push_back(std::cos(harmonic * phiTrack1) * pTTrack1 / ampl);
    tracksQy.push_back(std::sin(harmonic * phiTrack1) * pTTrack1 / ampl);
    tracksQx.push_back(std::cos(harmonic * phiTrack2) * pTTrack2 / ampl);
    tracksQy.push_back(std::sin(harmonic * phiTrack2) * pTTrack2 / ampl);
    tracksQx.push_back(std::cos(harmonic * phiTrack3) * pTTrack3 / ampl);
    tracksQy.push_back(std::sin(harmonic * phiTrack3) * pTTrack3 / ampl);
  }
  /// Compute the delta psi in the range [0, pi/harmonic]
  /// \param psi1 is the first angle
  /// \param psi2 is the second angle
  /// \note Ported from AliAnalysisTaskSECharmHadronvn::GetDeltaPsiSubInRange
  float getDeltaPsiInRange(float psi1, float psi2)
  {
    float deltaPsi = psi1 - psi2;
    deltaPsi = RecoDecay::constrainAngle(deltaPsi, -o2::constants::math::PI / harmonic, harmonic);
    return deltaPsi;
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
               const float cent,
               const float cosNPhi,
               const float sinNPhi,
               const float cosDeltaPhi,
               const float sp,
               const std::vector<float>& outputMl,
               const float occupancy,
               const o2::hf_evsel::HfCollisionRejectionMask hfevselflag)
  {
    if (occEstimator != 0) {
      std::vector<int> evtSelFlags = getEventSelectionFlags(hfevselflag);
      if (storeMl) {
        if (storeEP) {
          registry.fill(HIST("hSparseFlowCharm"), mass, pt, cent, sp, cosNPhi, sinNPhi, cosDeltaPhi, outputMl[0], outputMl[1], occupancy,
                        evtSelFlags[0], evtSelFlags[1], evtSelFlags[2], evtSelFlags[3], evtSelFlags[4]);
        } else {
          registry.fill(HIST("hSparseFlowCharm"), mass, pt, cent, sp, outputMl[0], outputMl[1], occupancy,
                        evtSelFlags[0], evtSelFlags[1], evtSelFlags[2], evtSelFlags[3], evtSelFlags[4]);
        }
      } else {
        if (storeEP) {
          registry.fill(HIST("hSparseFlowCharm"), mass, pt, cent, sp, cosNPhi, sinNPhi, cosDeltaPhi, occupancy,
                        evtSelFlags[0], evtSelFlags[1], evtSelFlags[2], evtSelFlags[3], evtSelFlags[4]);
        } else {
          registry.fill(HIST("hSparseFlowCharm"), mass, pt, cent, sp, occupancy,
                        evtSelFlags[0], evtSelFlags[1], evtSelFlags[2], evtSelFlags[3], evtSelFlags[4]);
        }
      }
    } else {
      if (storeMl) {
        if (storeEP) {
          registry.fill(HIST("hSparseFlowCharm"), mass, pt, cent, sp, cosNPhi, sinNPhi, cosDeltaPhi, outputMl[0], outputMl[1]);
        } else {
          registry.fill(HIST("hSparseFlowCharm"), mass, pt, cent, sp, outputMl[0], outputMl[1]);
        }
      } else {
        if (storeEP) {
          registry.fill(HIST("hSparseFlowCharm"), mass, pt, cent, sp, cosNPhi, sinNPhi, cosDeltaPhi);
        } else {
          registry.fill(HIST("hSparseFlowCharm"), mass, pt, cent, sp);
        }
      }
    }
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

  /// Get the Q vector
  /// \param collision is the collision with the Q vector information
  std::vector<float> getQvec(CollsWithQvecs::iterator const& collision)
  {
    float xQVec = -999.;
    float yQVec = -999.;
    float amplQVec = -999.;
    switch (qvecDetector) {
      case QvecEstimator::FV0A:
        xQVec = collision.qvecFV0ARe();
        yQVec = collision.qvecFV0AIm();
        break;
      case QvecEstimator::FT0M:
        xQVec = collision.qvecFT0MRe();
        yQVec = collision.qvecFT0MIm();
        break;
      case QvecEstimator::FT0A:
        xQVec = collision.qvecFT0ARe();
        yQVec = collision.qvecFT0AIm();
        break;
      case QvecEstimator::FT0C:
        xQVec = collision.qvecFT0CRe();
        yQVec = collision.qvecFT0CIm();
        break;
      case QvecEstimator::TPCPos:
        xQVec = collision.qvecBPosRe();
        yQVec = collision.qvecBPosIm();
        amplQVec = collision.nTrkBPos();
        break;
      case QvecEstimator::TPCNeg:
        xQVec = collision.qvecBNegRe();
        yQVec = collision.qvecBNegIm();
        amplQVec = collision.nTrkBNeg();
        break;
      case QvecEstimator::TPCTot:
        xQVec = collision.qvecBTotRe();
        yQVec = collision.qvecBTotIm();
        amplQVec = collision.nTrkBTot();
        break;
      default:
        LOG(warning) << "Q vector estimator not valid. Please choose between FV0A, FT0M, FT0A, FT0C, TPC Pos, TPC Neg. Fallback to FV0A";
        xQVec = collision.qvecFV0ARe();
        yQVec = collision.qvecFV0AIm();
        break;
    }
    return {xQVec, yQVec, amplQVec};
  }

  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param candidates are the selected candidates
  template <DecayChannel Channel, typename T1>
  void runFlowAnalysis(CollsWithQvecs::iterator const& collision,
                       T1 const& candidates)
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

    std::vector<float> qVecs = getQvec(collision);
    float xQVec = qVecs[0];
    float yQVec = qVecs[1];
    float const amplQVec = qVecs[2];
    float const evtPl = epHelper.GetEventPlane(xQVec, yQVec, harmonic);
    int nProngs = 3;

    for (const auto& candidate : candidates) {
      float massCand = 0.;
      std::vector<float> outputMl = {-999., -999.};

      if constexpr (std::is_same_v<T1, CandDsData> || std::is_same_v<T1, CandDsDataWMl>) {
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
        if constexpr (std::is_same_v<T1, CandDplusDataWMl>) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbDplusToPiKPi()[classMl->at(iclass)];
          }
        }
      } else if constexpr (std::is_same_v<T1, CandD0Data> || std::is_same_v<T1, CandD0DataWMl>) {
        nProngs = 2;
        switch (Channel) {
          case DecayChannel::D0ToPiK:
            massCand = HfHelper::invMassD0ToPiK(candidate);
            if constexpr (std::is_same_v<T1, CandD0DataWMl>) {
              for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
                outputMl[iclass] = candidate.mlProbD0()[classMl->at(iclass)];
              }
            }
            break;
          case DecayChannel::D0ToKPi:
            massCand = HfHelper::invMassD0barToKPi(candidate);
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
        if constexpr (std::is_same_v<T1, CandXic0DataWMl>) {
          for (unsigned int iclass = 0; iclass < classMl->size(); iclass++) {
            outputMl[iclass] = candidate.mlProbToXiPi()[classMl->at(iclass)];
          }
        }
      }

      float ptCand = 0.;
      float phiCand = 0.;

      if constexpr (std::is_same_v<T1, CandXic0Data> || std::is_same_v<T1, CandXic0DataWMl>) {
        ptCand = RecoDecay::sqrtSumOfSquares(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
        phiCand = std::atan2(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());
      } else {
        ptCand = candidate.pt();
        phiCand = candidate.phi();
      }

      // If TPC is used for the SP estimation, the tracks of the hadron candidate must be removed from the TPC Q vector to avoid double counting
      if (qvecDetector == QvecEstimator::TPCNeg || qvecDetector == QvecEstimator::TPCPos) {
        float const ampl = amplQVec - static_cast<float>(nProngs);
        std::vector<float> tracksQx = {};
        std::vector<float> tracksQy = {};
        if constexpr (std::is_same_v<T1, CandXic0Data> || std::is_same_v<T1, CandXic0DataWMl>) {
          // std::cout<<candidate.pxProng0()<<std::endl;
          getQvecXic0Tracks(candidate, tracksQx, tracksQy, ampl);
        } else {
          getQvecDtracks<Channel>(candidate, tracksQx, tracksQy, ampl);
        }
        for (auto iTrack{0u}; iTrack < tracksQx.size(); ++iTrack) {
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
        fillThn(massCand, ptCand, cent, cosNPhi, sinNPhi, cosDeltaPhi, scalprodCand, outputMl, occupancy, hfevflag);
      }
    }
  }

  // Ds with ML
  void processDsMl(CollsWithQvecs::iterator const& collision,
                   CandDsDataWMl const&)
  {
    auto candsDsToKKPiWMl = selectedDsToKKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsDsToPiKKWMl = selectedDsToPiKKWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::DsToKKPi>(collision, candsDsToKKPiWMl);
    runFlowAnalysis<DecayChannel::DsToPiKK>(collision, candsDsToPiKKWMl);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDsMl, "Process Ds candidates with ML", false);

  // Ds with rectangular cuts
  void processDs(CollsWithQvecs::iterator const& collision,
                 CandDsData const&)
  {
    auto candsDsToKKPi = selectedDsToKKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsDsToPiKK = selectedDsToPiKK->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::DsToKKPi>(collision, candsDsToKKPi);
    runFlowAnalysis<DecayChannel::DsToPiKK>(collision, candsDsToPiKK);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDs, "Process Ds candidates", false);

  // Dplus with ML
  void processDplusMl(CollsWithQvecs::iterator const& collision,
                      CandDplusDataWMl const& candidatesDplus)
  {
    runFlowAnalysis<DecayChannel::DplusToPiKPi>(collision, candidatesDplus);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDplusMl, "Process Dplus candidates with ML", false);

  // Dplus with rectangular cuts
  void processDplus(CollsWithQvecs::iterator const& collision,
                    CandDplusData const& candidatesDplus)
  {
    runFlowAnalysis<DecayChannel::DplusToPiKPi>(collision, candidatesDplus);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processDplus, "Process Dplus candidates", true);

  // D0 with ML
  void processD0Ml(CollsWithQvecs::iterator const& collision,
                   CandD0DataWMl const&)
  {
    auto candsD0ToPiKWMl = selectedD0ToPiKWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsD0ToKPiWMl = selectedD0ToKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::D0ToPiK>(collision, candsD0ToPiKWMl);
    runFlowAnalysis<DecayChannel::D0ToKPi>(collision, candsD0ToKPiWMl);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processD0Ml, "Process D0 candidates with ML", false);

  // D0 with rectangular cuts
  void processD0(CollsWithQvecs::iterator const& collision,
                 CandD0Data const&)
  {
    auto candsD0ToPiK = selectedD0ToPiK->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsD0ToKPi = selectedD0ToKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::D0ToPiK>(collision, candsD0ToPiK);
    runFlowAnalysis<DecayChannel::D0ToKPi>(collision, candsD0ToKPi);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processD0, "Process D0 candidates", false);

  // Lc with ML
  void processLcMl(CollsWithQvecs::iterator const& collision,
                   CandLcDataWMl const&)
  {
    auto candsLcToPKPiWMl = selectedLcToPKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsLcToPiKPWMl = selectedLcToPiKPWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::LcToPKPi>(collision, candsLcToPKPiWMl);
    runFlowAnalysis<DecayChannel::LcToPiKP>(collision, candsLcToPiKPWMl);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processLcMl, "Process Lc candidates with ML", false);

  // Lc with rectangular cuts
  void processLc(CollsWithQvecs::iterator const& collision,
                 CandLcData const&)
  {
    auto candsLcToPKPi = selectedLcToPKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsLcToPiKP = selectedLcToPiKP->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::LcToPKPi>(collision, candsLcToPKPi);
    runFlowAnalysis<DecayChannel::LcToPiKP>(collision, candsLcToPiKP);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processLc, "Process Lc candidates", false);

  // Xic with ML
  void processXicMl(CollsWithQvecs::iterator const& collision,
                    CandXicDataWMl const&)
  {
    auto candsXicToPKPiWMl = selectedXicToPKPiWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsXicToPiKPWMl = selectedXicToPiKPWMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::XicToPKPi>(collision, candsXicToPKPiWMl);
    runFlowAnalysis<DecayChannel::XicToPiKP>(collision, candsXicToPiKPWMl);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processXicMl, "Process Xic candidates with ML", false);

  // Xic with rectangular cuts
  void processXic(CollsWithQvecs::iterator const& collision,
                  CandXicData const&)
  {
    auto candsXicToPKPi = selectedXicToPKPi->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    auto candsXicToPiKP = selectedXicToPiKP->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::XicToPKPi>(collision, candsXicToPKPi);
    runFlowAnalysis<DecayChannel::XicToPiKP>(collision, candsXicToPiKP);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processXic, "Process Xic candidates", false);

  // Xic0 with ML
  void processXic0Ml(CollsWithQvecs::iterator const& collision,
                     CandXic0DataWMl const&)
  {
    auto candsXic0WMl = selectedXic0WMl->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::Xic0ToXiPi>(collision, candsXic0WMl);
  }
  PROCESS_SWITCH(HfTaskFlowCharmHadrons, processXic0Ml, "Process Xic0 candidates with ML", false);

  // Xic0
  void processXic0(CollsWithQvecs::iterator const& collision,
                   CandXic0Data const&)
  {
    auto candsXic0 = selectedXic0->sliceByCached(aod::hf_cand::collisionId, collision.globalIndex(), cache);
    runFlowAnalysis<DecayChannel::Xic0ToXiPi>(collision, candsXic0);
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

      registry.fill(HIST("epReso/hEpResoFT0cFT0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epFT0a)));
      registry.fill(HIST("epReso/hEpResoFT0cFV0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epFV0a)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBTots)));
      registry.fill(HIST("epReso/hEpResoFT0aFV0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epFV0a)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBTots)));
      registry.fill(HIST("epReso/hEpResoFT0mFV0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epFV0a)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBTots)));
      registry.fill(HIST("epReso/hEpResoFV0aTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFV0aTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFV0aTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epBTots)));
      registry.fill(HIST("epReso/hEpResoTPCposTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epBPoss, epBNegs)));
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
