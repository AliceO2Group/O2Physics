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

///
/// \file centralityQa.cxx
/// \brief This task does dedicated centrality QA
/// \author ALICE Collaboration
///

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonDataFormat/BunchFilling.h>
#include <DataFormatsFIT/Triggers.h>
#include <DataFormatsParameters/GRPECSObject.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TFormula.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TString.h>

#include <array>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

using CentRun3 = soa::Join<aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs, aod::CentFT0CVariant1s, aod::CentFT0CVariant2s, aod::CentFT0MAnchorCols, aod::CentFT0MAnchorBCs, aod::CentFT0MOuterAs, aod::FITExtraMults, aod::CentNTPVs>;
using BCsWithMatching = soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>;

struct CentralityQa {
  HistogramRegistry histos{"histos"};
  std::map<std::string, HistPtr> histPointers;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  std::bitset<o2::constants::lhc::LHCMaxBunches> collidingBunch;

  int runNumber{};
  uint64_t startOfRunTimestamp{};
  TList* hCentralityObjects = nullptr;

  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {1000, 0, 1000}, "Multiplicity"};
  ConfigurableAxis axisMultiplicityPV{"axisMultiplicityPV", {1000, 0, 1000}, "Multiplicity PV"};
  ConfigurableAxis axisChannelAmplitude{"axisChannelAmplitude", {5000, 0, 5000}, "Channel Amplitude"};
  ConfigurableAxis axisCentrality{"axisCentrality", {101, 0.0f, 101.0f}, "Centrality (%)"};

  Configurable<bool> loopOverMcCollisionsForMcHist{"loopOverMcCollisionsForMcHist", false, "Fill MC histograms in a loop over MC collisions? If no, they will be filled in a loop over reconstructed collisions"};

  struct : ConfigurableGroup {
    std::string prefix = "eventSelections"; // JSON group name
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
    Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
    Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border (Run 3 only)"};
    Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border (Run 3 only)"};
    Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track (Run 3 only)"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference (Run 3 only)"};
    Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF (Run 3 only)"};
    Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD (Run 3 only)"};
    Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInROFStd{"requireNoCollInROFStd", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF with mult. above a certain threshold (Run 3 only)"};
    Configurable<bool> requireNoCollInROFStrict{"requireNoCollInROFStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF (Run 3 only)"};
    Configurable<bool> requireINEL0{"requireINEL0", true, "require INEL>0 event selection"};
    Configurable<bool> requireINEL1{"requireINEL1", false, "require INEL>1 event selection"};

    Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position"};

    Configurable<bool> useFT0CbasedOccupancy{"useFT0CbasedOccupancy", false, "Use sum of FT0-C amplitudes for estimating occupancy? (if not, use track-based definition)"};
    // fast check on occupancy
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};

    Configurable<bool> requireIsBBT0A{"requireIsBBT0A", false, "Require beam-beam collisions based on timing information in FT0A"};
    Configurable<bool> requireIsBBT0C{"requireIsBBT0C", false, "Require beam-beam collisions based on timing information in FT0C"};

    Configurable<bool> rejectMismatchedBCs{"rejectMismatchedBCs", false, "Reject collision with BC different from MC BC"};
    Configurable<bool> rejectMismatchedFoundBCs{"rejectMismatchedFoundBCs", false, "Reject collision with found BC different from MC BC"};

    // Run 2 specific event selections
    Configurable<bool> requireSel7{"requireSel7", true, "require sel7 event selection (Run 2 only: event selection decision based on V0A & V0C)"};
    Configurable<bool> requireINT7{"requireINT7", true, "require INT7 trigger selection (Run 2 only)"};
    Configurable<bool> rejectIncompleteDAQ{"rejectIncompleteDAQ", true, "reject events with incomplete DAQ (Run 2 only)"};
    Configurable<bool> requireConsistentSPDAndTrackVtx{"requireConsistentSPDAndTrackVtx", true, "reject events with inconsistent in SPD and Track vertices (Run 2 only)"};
    Configurable<bool> rejectPileupFromSPD{"rejectPileupFromSPD", true, "reject events with pileup according to SPD vertexer (Run 2 only)"};
    Configurable<bool> rejectV0PFPileup{"rejectV0PFPileup", false, "reject events tagged as OOB pileup according to V0 past-future info (Run 2 only)"};
    Configurable<bool> rejectPileupInMultBins{"rejectPileupInMultBins", true, "reject events tagged as pileup according to multiplicity-differential pileup checks (Run 2 only)"};
    Configurable<bool> rejectPileupMV{"rejectPileupMV", true, "reject events tagged as pileup according to according to multi-vertexer (Run 2 only)"};
    Configurable<bool> rejectTPCPileup{"rejectTPCPileup", false, "reject events tagged as pileup according to pileup in TPC (Run 2 only)"};
    Configurable<bool> requireNoV0MOnVsOffPileup{"requireNoV0MOnVsOffPileup", false, "reject events tagged as OOB pileup according to online-vs-offline VOM correlation (Run 2 only)"};
    Configurable<bool> requireNoSPDOnVsOffPileup{"requireNoSPDOnVsOffPileup", false, "reject events tagged as pileup according to online-vs-offline SPD correlation (Run 2 only)"};
    Configurable<bool> requireNoSPDClsVsTklBG{"requireNoSPDClsVsTklBG", true, "reject events tagged as beam-gas and pileup according to cluster-vs-tracklet correlation (Run 2 only)"};
  } eventSelections;

  struct : ConfigurableGroup {
    std::string prefix = "bcsel";
    Configurable<bool> selectCollidingBCs{"selectCollidingBCs", true, "select colliding BCs"};
    Configurable<bool> selectTVX{"selectTVX", true, "select TVX"};
    Configurable<bool> selectFV0OrA{"selectFV0OrA", true, "select FV0 or A"};
    Configurable<bool> selectVertexZwithT0{"selectVertexZwithT0", false, "select vertex Z with T0"};
    Configurable<float> vertexZwithT0{"vertexZwithT0", 1000, "vertex Z with T0"};
    Configurable<bool> selectBBT0{"selectBBT0", false, "select BBT0"};
    Configurable<bool> rejectZNAC{"rejectZNAC", false, "reject ZNAC"};
    Configurable<bool> rejectIsFlangeEvent{"rejectIsFlangeEvent", false, "reject is flange event"};
  } bcsel;

  struct : ConfigurableGroup {
    std::string prefix = "studies";
    Configurable<bool> fv0a{"fv0a", true, "Enable centrality QA for fv0a"};
    Configurable<bool> ft0m{"ft0m", true, "Enable centrality QA for ft0m"};
    Configurable<bool> ft0mOuterA{"ft0mOuterA", true, "Enable centrality QA for ft0m"};
    Configurable<bool> ft0mAnchorCol{"ft0mAnchorCol", true, "Enable centrality QA for ft0m"};
    Configurable<bool> ft0mAnchorBc{"ft0mAnchorBc", true, "Enable centrality QA for ft0m"};
    Configurable<bool> ft0a{"ft0a", false, "Enable centrality QA for ft0a"};
    Configurable<bool> ft0c{"ft0c", true, "Enable centrality QA for ft0c"};
    Configurable<bool> ft0cVar1{"ft0cVar1", false, "Enable centrality QA for ft0cVar1"};
    Configurable<bool> ft0cVar2{"ft0cVar2", false, "Enable centrality QA for ft0cVar2"};
    Configurable<bool> fddm{"fddm", false, "Enable centrality QA for fddm"};
    Configurable<bool> ntpv{"ntpv", false, "Enable centrality QA for ntpv"};
    Configurable<bool> nGlo{"nGlo", true, "Enable centrality QA for nGlo"};
    Configurable<bool> mft{"mft", true, "Enable centrality QA for mft"};
  } studies;

  struct : ConfigurableGroup {
    std::string prefix = "centrality";
    Configurable<bool> useCustomCalibration{"useCustomCalibration", false, "override the centrality from the central calibration with a different calibration provided in pathCentrality"};
    Configurable<std::string> ccdbURL{"ccdbURL", "http://alice-ccdb.cern.ch", "ccdb url"};
    Configurable<std::string> pathCentrality{"pathCentrality", "Centrality/Estimators", "path to centrality calibration if useCustomCalibration is enabled"};
    Configurable<std::string> generator{"generator", "", "E.g. PYTHIA"};
  } centrality;

  enum EstimatorIndex { FV0A,
                        FT0M,
                        FT0MAnchorCol,
                        FT0MAnchorBC,
                        FT0MOuterA,
                        FT0A,
                        FT0C,
                        FT0CVar1,
                        FT0CVar2,
                        FDDM,
                        NTPV,
                        NGlobal,
                        MFT,
                        NEstimators };

  static constexpr int NSuperCalibPars = 6;
  static constexpr float CentralityNotFound = 105.f;

  struct Estimator {
    CentralityQa* outer = nullptr;
    bool doStudy = false;
    std::string estName;
    std::string histName;
    std::array<float, NSuperCalibPars> mcScalePars{};
    TH1* hCentrality = nullptr;
    TFormula* mcScale = nullptr;
    Estimator() = default;
    explicit Estimator(CentralityQa* o, bool b, std::string s0, std::string s1) : outer(o), doStudy(b), estName(std::move(s0)), histName(std::move(s1)) {}
    float getCentrality(const float mult, const float centTable)
    {
      if (outer->centrality.useCustomCalibration) {
        float lMult = mult;
        if (outer->doprocessRun3MonteCarlo && outer->hCentralityObjects != nullptr) {
          mcScale = dynamic_cast<TFormula*>(outer->hCentralityObjects->FindObject(TString::Format("%s-%s", outer->centrality.generator.value.c_str(), estName.c_str()).Data()));
          if (!mcScale) {
            return CentralityNotFound;
          }

          for (int ixpar = 0; ixpar < NSuperCalibPars; ++ixpar) {
            mcScalePars[ixpar] = mcScale->GetParameter(ixpar);
          }

          auto scaleMC = [](float x, const std::array<float, NSuperCalibPars>& pars) {
            float core = ((pars[0] + pars[1] * std::pow(x, pars[2])) - pars[3]) / pars[4];
            if (core < 0.0f) {
              return 0.0f; // this should be marked as low multiplicity and not mapped, core^pars[5] would be NaN
            }
            return std::pow(core, 1.0f / pars[5]);
          };

          lMult = scaleMC(mult, mcScalePars);
        }
        return hCentrality ? hCentrality->GetBinContent(hCentrality->FindBin(lMult)) : CentralityNotFound;
      }
      return centTable;
    }
    template <typename TCollision>
    void configure(const TCollision& col)
    {
      if (!outer->centrality.useCustomCalibration) {
        return;
      }

      if (!col.has_foundBC()) {
        return;
      }

      const auto& bc = col.template foundBC_as<BCsWithMatching>();
      if (bc.runNumber() != outer->runNumber) {
        outer->runNumber = bc.runNumber();
        LOGF(info, "Acquiring centrality calibration for run %i", outer->runNumber);
        outer->hCentralityObjects = outer->ccdb->getForRun<TList>(outer->centrality.pathCentrality, outer->runNumber);
        if (!outer->hCentralityObjects) {
          LOGF(info, "No centrality calibration list found for run %i", outer->runNumber);
        }
      }

      if (!outer->hCentralityObjects) {
        hCentrality = nullptr;
        return;
      }

      hCentrality = dynamic_cast<TH1*>(outer->hCentralityObjects->FindObject(Form("hCalibZeq%s", estName.c_str())));
      if (!hCentrality) {
        LOGF(debug, "Calibration missing for %s", estName.c_str());
      } else {
        LOGF(debug, "Calibration loaded for %s", estName.c_str());
      }
    }
  };

  std::vector<Estimator> estimators;
  Estimator initEstimator(const bool doStudy, const std::string& estName, const std::string& histName)
  {
    return Estimator(this, doStudy, estName, histName);
  }

  template <typename... Args>
  void insertHist(const std::string& name, const std::string& title, HistType type, const std::vector<AxisSpec>& axes)
  {
    histPointers[name] = histos.add(name.c_str(), title.c_str(), type, axes);
  }

  template <typename T>
  std::shared_ptr<T>& getHist(const std::string& name)
  {
    return std::get<std::shared_ptr<T>>(histPointers[name]);
  }

  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, CentRun3>> perMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentNGlobals>> perMcCollisionNGlobal = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MFTMults, aod::MultsExtra, aod::CentMFTs>> perMcCollisionMFT = aod::mccollisionlabel::mcCollisionId;

  void init(o2::framework::InitContext& /*initContext*/)
  {
    ccdb->setURL(centrality.ccdbURL);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    estimators.resize(NEstimators);
    estimators[FV0A] = initEstimator(studies.fv0a.value, "FV0", "FV0A");
    estimators[FT0M] = initEstimator(studies.ft0m.value, "FT0", "FT0M");
    estimators[FT0MAnchorCol] = initEstimator(studies.ft0mAnchorCol.value, "FT0MAnchorCol", "FT0MAnchorCol");
    estimators[FT0MAnchorBC] = initEstimator(studies.ft0mAnchorBc.value, "FT0MAnchorBc", "FT0MAnchorBc");
    estimators[FT0MOuterA] = initEstimator(studies.ft0mOuterA.value, "FT0MOuterA", "FT0MOuterA");
    estimators[FT0A] = initEstimator(studies.ft0a.value, "FT0A", "FT0A");
    estimators[FT0C] = initEstimator(studies.ft0c.value, "FT0C", "FT0C");
    estimators[FT0CVar1] = initEstimator(studies.ft0cVar1.value, "FT0CVariant1", "FT0CVariant1");
    estimators[FT0CVar2] = initEstimator(studies.ft0cVar2.value, "FT0CVariant2", "FT0CVariant2");
    estimators[FDDM] = initEstimator(studies.fddm.value, "FDDM", "FDDM");
    estimators[NTPV] = initEstimator(studies.ntpv.value, "NTracksPV", "NTPV");
    estimators[NGlobal] = initEstimator(studies.nGlo.value, "NGlobal", "NGlobal");
    estimators[MFT] = initEstimator(studies.mft.value, "MFT", "MFT");

    if (doprocessRun2) {
      insertHist("hCentRun2V0M", ";V0M centrality (%)", kTH1D, {axisCentrality});
      insertHist("hCentRun2V0A", ";V0A centrality (%)", kTH1D, {axisCentrality});
      insertHist("hCentRun2SPDTks", ";SPD tracklet centrality (%)", kTH1D, {axisCentrality});
      insertHist("hCentRun2SPDCls", ";SPD cluster centrality (%)", kTH1D, {axisCentrality});
      insertHist("hCentRun2CL0", ";CL0 centrality (%)", kTH1D, {axisCentrality});
      insertHist("hCentRun2CL1", ";CL1 centrality (%)", kTH1D, {axisCentrality});
    }

    if (doprocessRun3 || doprocessRun3MonteCarlo) {
      for (int iEst = 0; iEst < NEstimators; ++iEst) {
        const Estimator& est = estimators[iEst];
        if (!est.doStudy) {
          continue;
        }
        insertHist("hCent" + est.histName, ";" + est.histName + " centrality (%)", kTH1D, {axisCentrality});
        insertHist("hCentProfile" + est.histName, ";" + est.histName + " centrality (%)", kTProfile, {axisCentrality});
        insertHist("hMultEta05VsCent" + est.histName, ";" + est.histName + " Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisCentrality, axisMultiplicityPV});
        if (doprocessRun3MonteCarlo) {
          insertHist("hMultEta05VsGenMult" + est.histName, ";Multiplicity " + est.histName + ";Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
          insertHist("hGenMultEta05VsCentrality" + est.histName, ";" + est.histName + " Centrality (%); Generated multiplicity (|#it{#eta}| < 0.5)", kTH2D, {axisCentrality, axisMultiplicityPV});
          insertHist("hGenMultVsCentrality" + est.histName, ";" + est.histName + " Centrality (%); Generated multiplicity " + est.estName, kTH2D, {axisCentrality, axisMultiplicityPV});
        }
      }
    }

    if (doprocessBunchCrossings) {
      insertHist("hBCSelection", "hBCSelection", kTH1D, {{20, -0.5, 19.5f}});
      getHist<TH1>("hBCSelection")->GetXaxis()->SetBinLabel(1, "All BCs");
      getHist<TH1>("hBCSelection")->GetXaxis()->SetBinLabel(2, "Colliding BCs");
      getHist<TH1>("hBCSelection")->GetXaxis()->SetBinLabel(3, "TVX");
      getHist<TH1>("hBCSelection")->GetXaxis()->SetBinLabel(4, "FV0OrA");
      getHist<TH1>("hBCSelection")->GetXaxis()->SetBinLabel(5, "FT0PosZ");
      getHist<TH1>("hBCSelection")->GetXaxis()->SetBinLabel(6, "BB with FT0");
      getHist<TH1>("hBCSelection")->GetXaxis()->SetBinLabel(7, "zdc rej");
      getHist<TH1>("hBCSelection")->GetXaxis()->SetBinLabel(8, "isFlangeEvent");
      insertHist("hAmpVsChFT0A", "hAmpVsChFT0A;Channel; Amplitude", kTH2D, {{96, -0.5, 95.5}, axisChannelAmplitude});
      insertHist("hAmpVsChFT0C", "hAmpVsChFT0C;Channel; Amplitude", kTH2D, {{112, -0.5, 111.5}, axisChannelAmplitude});
    }

    histos.print();
  }

  template <typename TCollision>
  bool isCollisionAccepted(TCollision const& collision)
  // check whether the collision passes our collision selections
  {
    if (eventSelections.requireSel8 && !collision.sel8()) {
      return false;
    }

    if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }

    if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }

    if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }

    if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }

    if (eventSelections.requireIsBBT0A && !collision.selection_bit(aod::evsel::kIsBBT0A)) {
      return false;
    }

    if (eventSelections.requireIsBBT0C && !collision.selection_bit(aod::evsel::kIsBBT0C)) {
      return false;
    }

    if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }

    if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }

    if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }

    if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }

    if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }

    if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }

    if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }

    if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }

    if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }

    if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }

    static constexpr int OneTrackInEta1 = 1;
    if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < OneTrackInEta1) {
      return false;
    }

    static constexpr int TwoTracksInEta1 = 2;
    if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < TwoTracksInEta1) {
      return false;
    }

    float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
    if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
      return false;
    }

    if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
      return false;
    }

    if constexpr (requires { collision.has_mcCollision(); }) { // check if we are in MC
      if (!collision.has_mcCollision()) {
        return false;
      }

      const auto& mcCollision = collision.template mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
      const auto& recoBC = collision.template bc_as<BCsWithMatching>();
      const auto& foundBC = collision.template foundBC_as<BCsWithMatching>();
      const auto& mcBC = mcCollision.template bc_as<BCsWithMatching>();

      // Check that the BC in data and MC is the same
      if (eventSelections.rejectMismatchedBCs && recoBC.globalBC() != mcBC.globalBC()) {
        return false;
      }
      if (eventSelections.rejectMismatchedFoundBCs && foundBC.globalBC() != mcBC.globalBC()) {
        return false;
      }
    }

    return true;
  }

  template <typename TCollision>
  bool isCollisionAcceptedRun2(TCollision const& collision)
  // check whether the collision passes our collision selections
  {
    if (eventSelections.requireSel8 && !collision.sel8()) {
      return false;
    }

    if (eventSelections.requireSel7 && !collision.sel7()) {
      return false;
    }

    if (eventSelections.requireINT7 && !collision.alias_bit(kINT7)) {
      return false;
    }

    if (eventSelections.requireTriggerTVX && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }

    if (eventSelections.rejectIncompleteDAQ && !collision.selection_bit(o2::aod::evsel::kNoIncompleteDAQ)) {
      return false;
    }

    if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }

    if (eventSelections.requireConsistentSPDAndTrackVtx && !collision.selection_bit(o2::aod::evsel::kNoInconsistentVtx)) {
      return false;
    }

    if (eventSelections.rejectPileupFromSPD && !collision.selection_bit(o2::aod::evsel::kNoPileupFromSPD)) {
      return false;
    }

    if (eventSelections.rejectV0PFPileup && !collision.selection_bit(o2::aod::evsel::kNoV0PFPileup)) {
      return false;
    }

    if (eventSelections.rejectPileupInMultBins && !collision.selection_bit(o2::aod::evsel::kNoPileupInMultBins)) {
      return false;
    }

    if (eventSelections.rejectPileupMV && !collision.selection_bit(o2::aod::evsel::kNoPileupMV)) {
      return false;
    }

    if (eventSelections.rejectTPCPileup && !collision.selection_bit(o2::aod::evsel::kNoPileupTPC)) {
      return false;
    }

    if (eventSelections.requireNoV0MOnVsOffPileup && !collision.selection_bit(o2::aod::evsel::kNoV0MOnVsOfPileup)) {
      return false;
    }

    if (eventSelections.requireNoSPDOnVsOffPileup && !collision.selection_bit(o2::aod::evsel::kNoSPDOnVsOfPileup)) {
      return false;
    }

    if (eventSelections.requireNoSPDClsVsTklBG && !collision.selection_bit(o2::aod::evsel::kNoSPDClsVsTklBG)) {
      return false;
    }

    static constexpr int OneTrackInEta1 = 1;
    if (eventSelections.requireINEL0 && collision.multNTracksPVeta1() < OneTrackInEta1) {
      return false;
    }

    static constexpr int TwoTracksInEta1 = 2;
    if (eventSelections.requireINEL1 && collision.multNTracksPVeta1() < TwoTracksInEta1) {
      return false;
    }

    return true;
  }

  template <typename TBunchCrossing>
  bool isBunchCrossingAccepted(const TBunchCrossing& bc, bool fillHistograms = false)
  {
    if (fillHistograms) {
      histos.fill(HIST("hBCSelection"), 0); // all BCs
    }

    if (bc.runNumber() != runNumber) {
      runNumber = bc.runNumber();
      auto grpo = ccdb->getForRun<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", runNumber);
      startOfRunTimestamp = grpo->getTimeStart();
      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", startOfRunTimestamp);
      collidingBunch = grplhcif->getBunchFilling().getBCPattern();
    }

    const int localBC = bc.globalBC() % o2::constants::lhc::LHCMaxBunches;
    const bool collidingBC = collidingBunch.test(localBC);

    if (bcsel.selectCollidingBCs && !collidingBC) {
      return false;
    }

    if (fillHistograms) {
      histos.fill(HIST("hBCSelection"), 1); // colliding
    }

    if (bcsel.selectTVX && !bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }

    if (fillHistograms) {
      histos.fill(HIST("hBCSelection"), 2); // TVX
    }

    bool isFV0OrA = false;
    if (bc.has_fv0a()) {
      const auto& fv0 = bc.fv0a();
      std::bitset<8> fv0TriggerMask = fv0.triggerMask();
      isFV0OrA = fv0TriggerMask[o2::fit::Triggers::bitA];
    }

    if (bcsel.selectFV0OrA && !isFV0OrA) {
      return false;
    }

    if (fillHistograms) {
      histos.fill(HIST("hBCSelection"), 3); // FV0OrA
    }

    const float largeVertexZ = 100.0f;
    if (bcsel.selectVertexZwithT0 && bcsel.vertexZwithT0 < largeVertexZ) {
      if (bc.has_ft0()) {
        const auto& ft0 = bc.ft0();
        if (!ft0.isValidTime()) {
          return false;
        }
        if (std::abs(ft0.posZ()) > bcsel.vertexZwithT0) {
          return false;
        }
      }
    }

    if (fillHistograms) {
      histos.fill(HIST("hBCSelection"), 4); // FT0PosZ
    }

    if (bcsel.selectBBT0 && !bc.selection_bit(o2::aod::evsel::kIsBBT0A) && !bc.selection_bit(o2::aod::evsel::kIsBBT0C)) {
      return false;
    }

    if (fillHistograms) {
      histos.fill(HIST("hBCSelection"), 5); // t0ac time
    }

    if (bcsel.rejectZNAC && !bc.selection_bit(o2::aod::evsel::kIsBBZNA) && !bc.selection_bit(o2::aod::evsel::kIsBBZNC)) {
      return false;
    }

    if (fillHistograms) {
      histos.fill(HIST("hBCSelection"), 6); // znac time
    }

    if (bcsel.rejectIsFlangeEvent) {
      if (bc.has_ft0()) {
        const auto& ft0 = bc.ft0();
        constexpr int IsFlangeEventId = 7;
        std::bitset<8> ft0TriggerMask = ft0.triggerMask();
        if (ft0TriggerMask[IsFlangeEventId]) {
          return false;
        }
      }
    }

    if (fillHistograms) {
      histos.fill(HIST("hBCSelection"), 7); // isFlangeEvent
    }

    return true;
  }

  template <typename TCollision>
  void fillEstimatorHistos(const TCollision& col, Estimator est, const float mult, const float refCent)
  {
    if (!est.doStudy) {
      return;
    }

    est.configure(col);
    const auto cent = est.getCentrality(mult, refCent);
    LOGF(debug, "cent%s=%.0f", est.histName.c_str(), cent);
    getHist<TH1>("hCent" + est.histName)->Fill(cent);
    getHist<TH2>("hMultEta05VsCent" + est.histName)->Fill(cent, col.multNTracksPVetaHalf());
    getHist<TProfile>("hCentProfile" + est.histName)->Fill(cent, col.multNTracksPVetaHalf());
  }

  template <typename TMcCollision>
  void fillEstimatorMonteCarloHistos(const TMcCollision& mcCol, const Estimator& est, const float multMC, const float cent, const float nPVContribsInEta05)
  {
    if (!est.doStudy) {
      return;
    }

    getHist<TH2>("hMultEta05VsGenMult" + est.histName)->Fill(multMC, nPVContribsInEta05);
    getHist<TH2>("hGenMultEta05VsCentrality" + est.histName)->Fill(cent, mcCol.multMCNParticlesEta05());
    getHist<TH2>("hGenMultVsCentrality" + est.histName)->Fill(cent, multMC);
  }

  void processRun2(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2V0As, aod::CentRun2SPDTrks, aod::CentRun2SPDClss, aod::CentRun2CL0s, aod::CentRun2CL1s, aod::Mults>::iterator const& col)
  {
    if (!isCollisionAcceptedRun2(col)) {
      return;
    }
    LOGF(debug, "centV0M=%.0f", col.centRun2V0M());
    LOGF(debug, "centSPDTracklets=%.0f", col.centRun2SPDTracklets());
    LOGF(debug, "centSPDClusters=%.0f", col.centRun2SPDClusters());
    LOGF(debug, "centCL0=%.0f", col.centRun2CL0());
    LOGF(debug, "centCL1=%.0f", col.centRun2CL1());
    LOGF(debug, "centV0A=%.0f", col.centRun2V0A());

    getHist<TH1>("hCentRun2V0M")->Fill(col.centRun2V0M());
    getHist<TH1>("hCentRun2SPDTks")->Fill(col.centRun2SPDTracklets());
    getHist<TH1>("hCentRun2SPDCls")->Fill(col.centRun2SPDClusters());
    getHist<TH1>("hCentRun2CL0")->Fill(col.centRun2CL0());
    getHist<TH1>("hCentRun2CL1")->Fill(col.centRun2CL1());
    getHist<TH1>("hCentRun2V0A")->Fill(col.centRun2V0A());
  }

  void processRun3(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, CentRun3>::iterator const& col,
                   soa::Join<aod::MultsGlobal, aod::CentNGlobals> const& nGlobals,
                   soa::Join<aod::MFTMults, aod::CentMFTs> const& nMFTs,
                   BCsWithMatching const&)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    fillEstimatorHistos(col, estimators[FV0A], col.multFV0A(), col.centFV0A());
    fillEstimatorHistos(col, estimators[FT0M], col.multFT0M(), col.centFT0M());
    fillEstimatorHistos(col, estimators[FT0MAnchorCol], col.multFT0M(), col.centFT0MAnchorCol());
    fillEstimatorHistos(col, estimators[FT0MAnchorBC], col.multFT0M(), col.centFT0MAnchorBC());
    fillEstimatorHistos(col, estimators[FT0MOuterA], col.multFT0AOuter() + col.multFT0C(), col.centFT0MOuterA());
    fillEstimatorHistos(col, estimators[FT0A], col.multFT0A(), col.centFT0A());
    fillEstimatorHistos(col, estimators[FT0C], col.multFT0C(), col.centFT0C());
    fillEstimatorHistos(col, estimators[FT0CVar1], col.multFT0C(), col.centFT0CVariant1());
    fillEstimatorHistos(col, estimators[FT0CVar2], col.multFT0C(), col.centFT0CVariant2());
    fillEstimatorHistos(col, estimators[FDDM], col.multFDDM(), col.centFDDM());
    fillEstimatorHistos(col, estimators[NTPV], col.multNTracksPV(), col.centNTPV());

    if (nGlobals.size() > 0) {
      const auto& multCentNGlo = nGlobals.rawIteratorAt(col.globalIndex());
      fillEstimatorHistos(col, estimators[NGlobal], multCentNGlo.multNTracksGlobal(), multCentNGlo.centNGlobal());
    }

    if (nMFTs.size() > 0) {
      const auto& multCentNMFT = nMFTs.rawIteratorAt(col.globalIndex());
      fillEstimatorHistos(col, estimators[MFT], multCentNMFT.mftNtracks(), multCentNMFT.centMFT());
    }
  }

  void processRun3MonteCarlo(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, CentRun3> const& collisions,
                             soa::Join<aod::MultsGlobal, aod::CentNGlobals> const& nGlobals,
                             soa::Join<aod::MFTMults, aod::CentMFTs> const& nMFTs,
                             soa::Join<aod::McCollisions, aod::MultMCExtras> const& mcCollisions,
                             BCsWithMatching const& /*bcs*/)
  {
    for (const auto& col : collisions) {
      if (!isCollisionAccepted(col)) {
        continue;
      }

      fillEstimatorHistos(col, estimators[FV0A], col.multFV0A(), col.centFV0A());
      fillEstimatorHistos(col, estimators[FT0M], col.multFT0M(), col.centFT0M());
      fillEstimatorHistos(col, estimators[FT0MAnchorCol], col.multFT0M(), col.centFT0MAnchorCol());
      fillEstimatorHistos(col, estimators[FT0MAnchorBC], col.multFT0M(), col.centFT0MAnchorBC());
      fillEstimatorHistos(col, estimators[FT0MOuterA], col.multFT0AOuter() + col.multFT0C(), col.centFT0MOuterA());
      fillEstimatorHistos(col, estimators[FT0A], col.multFT0A(), col.centFT0A());
      fillEstimatorHistos(col, estimators[FT0C], col.multFT0C(), col.centFT0C());
      fillEstimatorHistos(col, estimators[FT0CVar1], col.multFT0C(), col.centFT0CVariant1());
      fillEstimatorHistos(col, estimators[FT0CVar2], col.multFT0C(), col.centFT0CVariant2());
      fillEstimatorHistos(col, estimators[FDDM], col.multFDDM(), col.centFDDM());
      fillEstimatorHistos(col, estimators[NTPV], col.multNTracksPV(), col.centNTPV());

      if (nGlobals.size() > 0) {
        const auto& multCentNGlo = nGlobals.rawIteratorAt(col.globalIndex());
        fillEstimatorHistos(col, estimators[NGlobal], multCentNGlo.multNTracksGlobal(), multCentNGlo.centNGlobal());
      }

      if (nMFTs.size() > 0) {
        const auto& multCentNMFT = nMFTs.rawIteratorAt(col.globalIndex());
        fillEstimatorHistos(col, estimators[MFT], multCentNMFT.mftNtracks(), multCentNMFT.centMFT());
      }

      if (!loopOverMcCollisionsForMcHist) {
        const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
        std::vector<float> cent;
        cent.resize(NEstimators, CentralityNotFound);

        cent[FV0A] = estimators[FV0A].getCentrality(col.multFV0A(), col.centFV0A());
        cent[FT0M] = estimators[FT0M].getCentrality(col.multFT0M(), col.centFT0M());
        cent[FT0MAnchorCol] = estimators[FT0MAnchorCol].getCentrality(col.multFT0M(), col.centFT0MAnchorCol());
        cent[FT0MAnchorBC] = estimators[FT0MAnchorBC].getCentrality(col.multFT0M(), col.centFT0MAnchorBC());
        cent[FT0MOuterA] = estimators[FT0MOuterA].getCentrality(col.multFT0AOuter() + col.multFT0C(), col.centFT0MOuterA());
        cent[FT0A] = estimators[FT0A].getCentrality(col.multFT0A(), col.centFT0A());
        cent[FT0C] = estimators[FT0C].getCentrality(col.multFT0C(), col.centFT0C());
        cent[FT0CVar1] = estimators[FT0CVar1].getCentrality(col.multFT0C(), col.centFT0CVariant1());
        cent[FT0CVar2] = estimators[FT0CVar2].getCentrality(col.multFT0C(), col.centFT0CVariant2());
        cent[FDDM] = estimators[FDDM].getCentrality(col.multFDDM(), col.centFDDM());
        cent[NTPV] = estimators[NTPV].getCentrality(col.multNTracksPV(), col.centNTPV());

        fillEstimatorMonteCarloHistos(mcCol, estimators[FV0A], mcCol.multMCFV0A(), cent[FV0A], col.multNTracksPVetaHalf());
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0M], mcCol.multMCFT0A() + mcCol.multMCFT0C(), cent[FT0M], col.multNTracksPVetaHalf());
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0MAnchorCol], mcCol.multMCFT0A() + mcCol.multMCFT0C(), cent[FT0MAnchorCol], col.multNTracksPVetaHalf());
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0MAnchorBC], mcCol.multMCFT0A() + mcCol.multMCFT0C(), cent[FT0MAnchorBC], col.multNTracksPVetaHalf());
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0MOuterA], mcCol.multMCFT0A() + mcCol.multMCFT0C(), cent[FT0MOuterA], col.multNTracksPVetaHalf());
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0A], mcCol.multMCFT0A(), cent[FT0A], col.multNTracksPVetaHalf());
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0C], mcCol.multMCFT0C(), cent[FT0C], col.multNTracksPVetaHalf());
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0CVar1], mcCol.multMCFT0C(), cent[FT0CVar1], col.multNTracksPVetaHalf());
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0CVar2], mcCol.multMCFT0C(), cent[FT0CVar2], col.multNTracksPVetaHalf());
        fillEstimatorMonteCarloHistos(mcCol, estimators[FDDM], mcCol.multMCFDDA() + mcCol.multMCFDDC(), cent[FDDM], col.multNTracksPVetaHalf());
        fillEstimatorMonteCarloHistos(mcCol, estimators[NTPV], mcCol.multMCNParticlesEta08(), cent[NTPV], col.multNTracksPVetaHalf());
        fillEstimatorMonteCarloHistos(mcCol, estimators[NGlobal], mcCol.multMCNParticlesEta08(), cent[NGlobal], col.multNTracksPVetaHalf());

        if (nGlobals.size() > 0) {
          const auto& multCentNGlo = nGlobals.rawIteratorAt(col.globalIndex());
          cent[NGlobal] = estimators[NGlobal].getCentrality(multCentNGlo.multNTracksGlobal(), multCentNGlo.centNGlobal());
          fillEstimatorMonteCarloHistos(mcCol, estimators[NGlobal], mcCol.multMCNParticlesEta08(), cent[NGlobal], col.multNTracksPVetaHalf());
        }

        if (nMFTs.size() > 0) {
          const auto& multCentNMFT = nMFTs.rawIteratorAt(col.globalIndex());
          cent[MFT] = estimators[MFT].getCentrality(multCentNMFT.mftNtracks(), multCentNMFT.centMFT());
          // fillEstimatorMonteCarloHistos(mcCol, estimators[MFT], mcCol.multMCMFT(), cent[MFT], col.multNTracksPVetaHalf()); // FIXME: uncomment when MC MFT mult is added in aod::MultMCExtras
        }
      }
    }

    if (loopOverMcCollisionsForMcHist) {
      for (const auto& mcCol : mcCollisions) {
        auto groupedCollisions = collisions.sliceBy(perMcCollision, mcCol.globalIndex());
        auto groupedCollisionsNGlobal = collisions.sliceBy(perMcCollisionNGlobal, mcCol.globalIndex());
        auto groupedCollisionsMFT = collisions.sliceBy(perMcCollisionMFT, mcCol.globalIndex());

        // Check if there is at least one of the reconstructed collisions associated to this MC collision
        // If so, we consider it
        int biggestNContribs = -1;
        int nContribsInEta05 = -1;
        std::vector<float> cent;
        cent.resize(NEstimators, CentralityNotFound);

        for (auto const& col : groupedCollisions) {
          if (!isCollisionAccepted(col)) {
            continue;
          }

          for (int iEst = 0; iEst < NEstimators; ++iEst) {
            if (iEst == NGlobal || iEst == MFT) {
              continue;
            }
            estimators[iEst].configure(col);
          }

          const float recoColCentFV0A = estimators[FV0A].getCentrality(col.multFV0A(), col.centFV0A());
          const float recoColCentFT0M = estimators[FT0M].getCentrality(col.multFT0M(), col.centFT0M());
          const float recoColCentFT0MAnchorCol = estimators[FT0MAnchorCol].getCentrality(col.multFT0M(), col.centFT0MAnchorCol());
          const float recoColCentFT0MAnchorBC = estimators[FT0MAnchorBC].getCentrality(col.multFT0M(), col.centFT0MAnchorBC());
          const float recoColCentFT0MOuterA = estimators[FT0MOuterA].getCentrality(col.multFT0AOuter() + col.multFT0C(), col.centFT0MOuterA());
          const float recoColCentFT0A = estimators[FT0A].getCentrality(col.multFT0A(), col.centFT0A());
          const float recoColCentFT0C = estimators[FT0C].getCentrality(col.multFT0C(), col.centFT0C());
          const float recoColCentFT0CVar1 = estimators[FT0CVar1].getCentrality(col.multFT0C(), col.centFT0CVariant1());
          const float recoColCentFT0CVar2 = estimators[FT0CVar2].getCentrality(col.multFT0C(), col.centFT0CVariant2());
          const float recoColCentFDDM = estimators[FDDM].getCentrality(col.multFDDM(), col.centFDDM());
          const float recoColCentNTPV = estimators[NTPV].getCentrality(col.multNTracksPV(), col.centNTPV());

          // One McCollision can be reconstructed multiple times
          // Use centrality from reconstructed collision with most PV contributors
          if (biggestNContribs < col.multPVTotalContributors()) {
            biggestNContribs = col.multPVTotalContributors();
            nContribsInEta05 = col.multNTracksPVetaHalf();
            cent[FV0A] = recoColCentFV0A;
            cent[FT0M] = recoColCentFT0M;
            cent[FT0MAnchorCol] = recoColCentFT0MAnchorCol;
            cent[FT0MAnchorBC] = recoColCentFT0MAnchorBC;
            cent[FT0MOuterA] = recoColCentFT0MOuterA;
            cent[FT0A] = recoColCentFT0A;
            cent[FT0C] = recoColCentFT0C;
            cent[FT0CVar1] = recoColCentFT0CVar1;
            cent[FT0CVar2] = recoColCentFT0CVar2;
            cent[FDDM] = recoColCentFDDM;
            cent[NTPV] = recoColCentNTPV;
          }
        }

        if (nGlobals.size() > 0 && studies.nGlo.value) {
          for (auto const& col : groupedCollisionsNGlobal) {
            estimators[NGlobal].configure(col);
            const auto& multCentNGlo = nGlobals.rawIteratorAt(col.globalIndex());
            const float recoColCentNGlo = estimators[NGlobal].getCentrality(multCentNGlo.multNTracksGlobal(), multCentNGlo.centNGlobal());
            if (biggestNContribs < col.multPVTotalContributors()) {
              biggestNContribs = col.multPVTotalContributors();
              nContribsInEta05 = col.multNTracksPVetaHalf();
              cent[NGlobal] = recoColCentNGlo;
            }
          }
        }

        if (nMFTs.size() > 0 && studies.mft.value) {
          for (auto const& col : groupedCollisionsMFT) {
            estimators[MFT].configure(col);
            const auto& multCentNMFT = nMFTs.rawIteratorAt(col.globalIndex());
            const float recoColCentNMFT = estimators[MFT].getCentrality(multCentNMFT.mftNtracks(), multCentNMFT.centMFT());
            if (biggestNContribs < col.multPVTotalContributors()) {
              biggestNContribs = col.multPVTotalContributors();
              nContribsInEta05 = col.multNTracksPVetaHalf();
              cent[MFT] = recoColCentNMFT;
            }
          }
        }

        fillEstimatorMonteCarloHistos(mcCol, estimators[FV0A], mcCol.multMCFV0A(), cent[FV0A], nContribsInEta05);
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0M], mcCol.multMCFT0A() + mcCol.multMCFT0C(), cent[FT0M], nContribsInEta05);
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0MAnchorCol], mcCol.multMCFT0A() + mcCol.multMCFT0C(), cent[FT0MAnchorCol], nContribsInEta05);
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0MAnchorBC], mcCol.multMCFT0A() + mcCol.multMCFT0C(), cent[FT0MAnchorBC], nContribsInEta05);
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0MOuterA], mcCol.multMCFT0A() + mcCol.multMCFT0C(), cent[FT0MOuterA], nContribsInEta05);
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0A], mcCol.multMCFT0A(), cent[FT0A], nContribsInEta05);
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0C], mcCol.multMCFT0C(), cent[FT0C], nContribsInEta05);
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0CVar1], mcCol.multMCFT0C(), cent[FT0CVar1], nContribsInEta05);
        fillEstimatorMonteCarloHistos(mcCol, estimators[FT0CVar2], mcCol.multMCFT0C(), cent[FT0CVar2], nContribsInEta05);
        fillEstimatorMonteCarloHistos(mcCol, estimators[FDDM], mcCol.multMCFDDA() + mcCol.multMCFDDC(), cent[FDDM], nContribsInEta05);
        fillEstimatorMonteCarloHistos(mcCol, estimators[NTPV], mcCol.multMCNParticlesEta08(), cent[NTPV], nContribsInEta05);
        fillEstimatorMonteCarloHistos(mcCol, estimators[NGlobal], mcCol.multMCNParticlesEta08(), cent[NGlobal], nContribsInEta05);
        // fillEstimatorMonteCarloHistos(mcCol, estimators[MFT], mcCol.multMCMFT(), cent[MFT], nContribsInEta05); // FIXME: uncomment when MC MFT mult is added in aod::MultMCExtras
      }
    }
  }

  using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
  void processBunchCrossings(soa::Join<BCsWithRun3Matchings, aod::BCFlags, aod::BcSels>::iterator const& bc, aod::FT0s const&, aod::FV0As const&)
  {
    if (!isBunchCrossingAccepted(bc, true)) {
      return;
    }

    if (bc.has_ft0()) {
      const auto& ft0 = bc.ft0();
      for (size_t ii{0}; ii < ft0.channelA().size(); ++ii) {
        histos.fill(HIST("hAmpVsChFT0A"), ft0.channelA()[ii], ft0.amplitudeA()[ii]);
      }

      for (size_t ii{0}; ii < ft0.channelC().size(); ++ii) {
        histos.fill(HIST("hAmpVsChFT0C"), ft0.channelC()[ii], ft0.amplitudeC()[ii]);
      }
    }
  }

  PROCESS_SWITCH(CentralityQa, processRun2, "Process with Run 2 V0A", false);
  PROCESS_SWITCH(CentralityQa, processRun3, "Process with Run 3", true);
  PROCESS_SWITCH(CentralityQa, processRun3MonteCarlo, "Process with Run 3 with MC", false);
  PROCESS_SWITCH(CentralityQa, processBunchCrossings, "Process with Run 3 BC table", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityQa>(cfgc)};
}
