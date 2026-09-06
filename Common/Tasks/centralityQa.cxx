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
/// \author ALICE
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
#include <TString.h>

#include <array>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <utility>

using namespace o2;
using namespace o2::framework;

struct CentralityQa {
  HistogramRegistry histos{"histos"};
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  std::bitset<o2::constants::lhc::LHCMaxBunches> collidingBunch;

  bool isRun2 = false;
  bool isMC = false;
  int runNumber{};
  uint64_t startOfRunTimestamp{};
  TList* hCentralityObjects = nullptr;

  Configurable<int> nBins{"nBins", 1050, "number of bins"};
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
    std::string prefix = "centrality";
    Configurable<bool> useCustomCalibration{"useCustomCalibration", false, "override the centrality from the central calibration with a different calibration provided in pathCentrality"};
    Configurable<std::string> ccdbURL{"ccdbURL", "http://alice-ccdb.cern.ch", "ccdb url"};
    Configurable<std::string> pathCentrality{"pathCentrality", "Centrality/Estimators", "path to centrality calibration if useCustomCalibration is enabled"};
    Configurable<std::string> generator{"generator", "", "E.g. PYTHIA"};
  } centrality;

  static constexpr int NSuperCalibPars = 6;
  static constexpr float CentralityNotFound = 105.f;

  struct Estimator {
    CentralityQa* outer = nullptr;
    std::string name;
    std::array<float, NSuperCalibPars> mcScalePars{};
    TH1* hCentrality = nullptr;
    TFormula* mcScale = nullptr;
    explicit Estimator(CentralityQa* o, std::string s) : outer(o), name(std::move(s)) {}
    float getCentrality(const float mult, const float centTable)
    {
      if (outer->centrality.useCustomCalibration) {
        float lMult = mult;
        if (outer->isMC && outer->hCentralityObjects != nullptr) {
          mcScale = dynamic_cast<TFormula*>(outer->hCentralityObjects->FindObject(TString::Format("%s-%s", outer->centrality.generator.value.c_str(), name.c_str()).Data()));
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
  };

  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFV0As>> perMcCollisionFV0A = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFT0Ms>> perMcCollisionFT0M = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFT0As>> perMcCollisionFT0A = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFT0Cs>> perMcCollisionFT0C = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFT0CVariant1s>> perMcCollisionFT0CVar1 = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFT0CVariant2s>> perMcCollisionFT0CVar2 = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFDDMs>> perMcCollisionFDDM = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentNGlobals>> perMcCollisionNGlobal = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentNTPVs>> perMcCollisionNTPV = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MFTMults, aod::MultsExtra, aod::CentMFTs>> perMcCollisionMFT = aod::mccollisionlabel::mcCollisionId;

  void init(o2::framework::InitContext& /*initContext*/)
  {
    ccdb->setURL(centrality.ccdbURL);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    isRun2 = doprocessRun2PP ||
             doprocessRun2PPb ||
             doprocessRun2PbPb;

    isMC = doprocessMonteCarloRun3_FV0A ||
           doprocessMonteCarloRun3_FT0M ||
           doprocessMonteCarloRun3_FT0A ||
           doprocessMonteCarloRun3_FT0C ||
           doprocessMonteCarloRun3_FT0CVar1 ||
           doprocessMonteCarloRun3_FT0CVar2 ||
           doprocessMonteCarloRun3_MFT ||
           doprocessMonteCarloRun3_NGlobal ||
           doprocessMonteCarloRun3_NTPV ||
           doprocessMonteCarloRun3_FT0MAnchorCol ||
           doprocessMonteCarloRun3_FT0MAnchorBC;

    if (isRun2) {
      histos.add("hCentRun2V0M", ";V0M centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentRun2V0A", ";V0A centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentRun2SPDTks", ";SPD tracklet centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentRun2SPDCls", ";SPD cluster centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentRun2CL0", ";CL0 centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentRun2CL1", ";CL1 centrality (%)", kTH1D, {{nBins, 0, 105.}});
    } else {
      histos.add("hCentFV0A", ";FV0A centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0M", ";FT0M centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0MOuterA", ";FT0M centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0A", ";FT0A centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0C", ";FT0C centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0CVar1", ";FT0CVar1 centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0CVar2", ";FT0CVar2 centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFDDM", ";FDDM centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentNTPV", ";NTPV centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentNGlobal", ";NGlobal centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentMFT", ";MFT centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0MAnchorCols", ";FT0MAnchorCols centrality (%)", kTH1D, {{nBins, 0, 105.}});
      histos.add("hCentFT0MAnchorBCs", ";FT0MAnchorBCs centrality (%)", kTH1D, {{nBins, 0, 105.}});

      // profiles of midrapidity multiplicity density
      histos.add("hCentProfileFV0A", ";FV0A centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0M", ";FT0M centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0MOuterA", ";FT0M centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0A", ";FT0A centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0C", ";FT0C centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0CVar1", ";FT0CVar1 centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0CVar2", ";FT0CVar2 centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFDDM", ";FDDM centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileNTPV", ";NTPV centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileNGlobal", ";NGlobal centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileMFT", ";MFT centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0MAnchorCols", ";FT0MAnchorCols centrality (%)", kTProfile, {{nBins, 0, 105.}});
      histos.add("hCentProfileFT0MAnchorBCs", ";FT0MAnchorBCs centrality (%)", kTProfile, {{nBins, 0, 105.}});

      histos.add("hMultEta05VsCentFV0A", ";FV0A centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0M", ";FT0M centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0MOuterA", ";FT0M centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0A", ";FT0A centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0C", ";FT0C centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0CVar1", ";FT0CVar1 centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0CVar2", ";FT0CVar2 centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFDDM", ";FDDM centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentNTPV", ";NTPV centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentNGlobal", ";NGlobal centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentMFT", ";MFT centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0MAnchorCols", ";FT0MAnchorCols centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});
      histos.add("hMultEta05VsCentFT0MAnchorBCs", ";FT0MAnchorBCs centrality (%); Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {{nBins, 0, 105.}, axisMultiplicityPV});

      if (isMC) {
        histos.add("hMultEta05VsGenMultFV0A", ";Multiplicity FV0A; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFT0M", ";Multiplicity FT0M; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFT0MAnchorCols", ";Multiplicity FT0MAnchorCols; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFT0MAnchorBCs", ";Multiplicity FT0MAnchorBCs; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFT0A", ";Multiplicity FT0A; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFT0C", ";Multiplicity FT0C; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFT0CVar1", ";Multiplicity FT0CVar1; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFT0CVar2", ";Multiplicity FT0CVar2; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultFDDM", ";Multiplicity FDDM; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultNTPV", ";Multiplicity NTPV; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultNGlobal", ";Multiplicity NGlobal; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});
        histos.add("hMultEta05VsGenMultMFT", ";Multiplicity MFT; Multiplicity PV contributors (|#it{#eta}| < 0.5)", kTH2D, {axisMultiplicity, axisMultiplicityPV});

        histos.add("hGenMultEta05VsCentralityFV0A", ";FV0A Centrality (%); Generated multiplicity (|#it{#eta}| < 0.5)", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultEta05VsCentralityFT0M", ";FT0M Centrality (%); Generated multiplicity (|#it{#eta}| < 0.5)", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultEta05VsCentralityFT0A", ";FT0A Centrality (%); Generated multiplicity (|#it{#eta}| < 0.5)", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultEta05VsCentralityFT0C", ";FT0C Centrality (%); Generated multiplicity (|#it{#eta}| < 0.5)", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultEta05VsCentralityFT0CVar1", ";FT0CVar1 Centrality (%); Generated multiplicity (|#it{#eta}| < 0.5)", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultEta05VsCentralityFT0CVar2", ";FT0CVar2 Centrality (%); Generated multiplicity (|#it{#eta}| < 0.5)", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultEta05VsCentralityFDDM", ";FDDM Centrality (%); Generated multiplicity (|#it{#eta}| < 0.5)", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultEta05VsCentralityNTPV", ";NTPV Centrality (%); Generated multiplicity (|#it{#eta}| < 0.5)", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultEta05VsCentralityNGlobal", ";NGlobal Centrality (%); Generated multiplicity (|#it{#eta}| < 0.5)", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultEta05VsCentralityMFT", ";MFT Centrality (%); Generated multiplicity (|#it{#eta}| < 0.5)", kTH2D, {axisCentrality, axisMultiplicityPV});

        histos.add("hGenMultVsCentralityFV0A", ";FV0A Centrality (%); Generated multiplicity FV0A", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultVsCentralityFT0M", ";FT0M Centrality (%); Generated multiplicity FT0M", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultVsCentralityFT0A", ";FT0A Centrality (%); Generated multiplicity FT0A", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultVsCentralityFT0C", ";FT0C Centrality (%); Generated multiplicity FT0C", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultVsCentralityFT0CVar1", ";FT0CVar1 Centrality (%); Generated multiplicity FT0C", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultVsCentralityFT0CVar2", ";FT0CVar2 Centrality (%); Generated multiplicity FT0C", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultVsCentralityFDDM", ";FDDM Centrality (%); Generated multiplicity FDDM", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultVsCentralityNTPV", ";NTPV Centrality (%); Generated multiplicity NTPV", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultVsCentralityNGlobal", ";NGlobal Centrality (%); Generated multiplicity global tracks", kTH2D, {axisCentrality, axisMultiplicityPV});
        histos.add("hGenMultVsCentralityMFT", ";MFT Centrality (%); Generated multiplicity MFT", kTH2D, {axisCentrality, axisMultiplicityPV});
      }
    }

    if (doprocessBunchCrossings) {
      histos.add("hBCSelection", "hBCSelection", kTH1D, {{20, -0.5, 19.5f}});
      histos.get<TH1>(HIST("hBCSelection"))->GetXaxis()->SetBinLabel(1, "All BCs");
      histos.get<TH1>(HIST("hBCSelection"))->GetXaxis()->SetBinLabel(2, "Colliding BCs");
      histos.get<TH1>(HIST("hBCSelection"))->GetXaxis()->SetBinLabel(3, "TVX");
      histos.get<TH1>(HIST("hBCSelection"))->GetXaxis()->SetBinLabel(4, "FV0OrA");
      histos.get<TH1>(HIST("hBCSelection"))->GetXaxis()->SetBinLabel(5, "FT0PosZ");
      histos.get<TH1>(HIST("hBCSelection"))->GetXaxis()->SetBinLabel(6, "BB with FT0");
      histos.get<TH1>(HIST("hBCSelection"))->GetXaxis()->SetBinLabel(7, "zdc rej");
      histos.get<TH1>(HIST("hBCSelection"))->GetXaxis()->SetBinLabel(8, "isFlangeEvent");
      histos.add("hAmpVsChFT0A", "hAmpVsChFT0A;Channel; Amplitude", kTH2D, {{96, -0.5, 95.5}, axisChannelAmplitude});
      histos.add("hAmpVsChFT0C", "hAmpVsChFT0C;Channel; Amplitude", kTH2D, {{112, -0.5, 111.5}, axisChannelAmplitude});
    }

    histos.print();
  }

  template <typename TCollision>
  Estimator initEstimator(const TCollision& col, const std::string& name)
  {
    Estimator est(this, name);
    if (centrality.useCustomCalibration) {
      if (!col.has_foundBC()) {
        return est;
      }

      const auto bc = col.template foundBC_as<aod::BCs>();
      if (bc.runNumber() != runNumber) {
        runNumber = bc.runNumber();
        LOGF(info, "Acquiring centrality calibration for run %i", runNumber);
        hCentralityObjects = ccdb->getForRun<TList>(centrality.pathCentrality, runNumber);
        if (!hCentralityObjects) {
          LOGF(info, "No centrality calibration list found for run %i", runNumber);
        }
      }

      if (!hCentralityObjects) {
        est.hCentrality = nullptr;
        return est;
      }

      est.hCentrality = dynamic_cast<TH1*>(hCentralityObjects->FindObject(Form("hCalibZeq%s", est.name.c_str())));
      if (!est.hCentrality) {
        LOGF(debug, "Calibration missing for %s", est.name.c_str());
      } else {
        LOGF(debug, "Calibration loaded for %s", est.name.c_str());
      }
    }

    return est;
  }

  template <typename TCollision>
  bool isCollisionAccepted(TCollision const& collision)
  // check whether the collision passes our collision selections
  {
    if constexpr (
      requires { collision.centFV0A(); } ||
      requires { collision.centFT0M(); } ||
      requires { collision.centFT0A(); } ||
      requires { collision.centFT0C(); } ||
      requires { collision.centFT0CVariant1(); } ||
      requires { collision.centFT0CVariant2(); } ||
      requires { collision.centFT0MOuterA(); } ||
      requires { collision.centFDDM(); } ||
      requires { collision.centNTPV(); } ||
      requires { collision.centNGlobal(); } ||
      requires { collision.centMFT(); } ||
      requires { collision.centFT0MAnchorCol(); } ||
      requires { collision.centFT0MAnchorBC(); }) { // check if we are in Run 3
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
        const auto& recoBC = collision.template bc_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>();
        const auto& foundBC = collision.template foundBC_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>();
        const auto& mcBC = mcCollision.template bc_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>();

        // Check that the BC in data and MC is the same
        if (eventSelections.rejectMismatchedBCs && recoBC.globalBC() != mcBC.globalBC()) {
          return false;
        }
        if (eventSelections.rejectMismatchedFoundBCs && foundBC.globalBC() != mcBC.globalBC()) {
          return false;
        }
      }
    } else { // we are in Run 2
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

  void processRun2PP(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2SPDTrks, aod::CentRun2SPDClss, aod::Mults>::iterator const& col)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }
    LOGF(debug, "centV0M=%.0f", col.centRun2V0M());
    LOGF(debug, "centSPDTracklets=%.0f", col.centRun2SPDTracklets());
    LOGF(debug, "centSPDClusters=%.0f", col.centRun2SPDClusters());
    // fill centrality histos
    histos.fill(HIST("hCentRun2V0M"), col.centRun2V0M());
    histos.fill(HIST("hCentRun2SPDTks"), col.centRun2SPDTracklets());
    histos.fill(HIST("hCentRun2SPDCls"), col.centRun2SPDClusters());
  }
  PROCESS_SWITCH(CentralityQa, processRun2PP, "Process with Run2 SPD clusters centrality/multiplicity estimation", false);

  void processRun2PbPb(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::CentRun2SPDTrks, aod::CentRun2CL0s, aod::CentRun2CL1s, aod::Mults>::iterator const& col)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }
    LOGF(debug, "centV0M=%.0f", col.centRun2V0M());
    LOGF(debug, "centSPDTracklets=%.0f", col.centRun2SPDTracklets());
    LOGF(debug, "centCL0=%.0f", col.centRun2CL0());
    LOGF(debug, "centCL1=%.0f", col.centRun2CL1());
    // fill centrality histos
    histos.fill(HIST("hCentRun2V0M"), col.centRun2V0M());
    histos.fill(HIST("hCentRun2SPDTks"), col.centRun2SPDTracklets());
    histos.fill(HIST("hCentRun2CL0"), col.centRun2CL0());
    histos.fill(HIST("hCentRun2CL1"), col.centRun2CL1());
  }
  PROCESS_SWITCH(CentralityQa, processRun2PbPb, "Process with Run2 CL0 and CL1 multiplicities centrality/multiplicity  estimation", false);

  void processRun2PPb(soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0As, aod::Mults>::iterator const& col)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }
    LOGF(debug, "centV0A=%.0f", col.centRun2V0A());
    // fill centrality histos
    histos.fill(HIST("hCentRun2V0A"), col.centRun2V0A());
  }
  PROCESS_SWITCH(CentralityQa, processRun2PPb, "Process with Run2 V0A multiplicitY centrality/multiplicity  estimation", false);

  void processRun3_FV0A(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::CentFV0As>::iterator const& col, aod::BCs const&)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator fv0a = initEstimator(col, "FV0");
    const float centFV0A = fv0a.getCentrality(col.multFV0A(), col.centFV0A());

    LOGF(debug, "centFV0A=%.0f", centFV0A);
    histos.fill(HIST("hCentFV0A"), centFV0A);
    histos.fill(HIST("hCentProfileFV0A"), centFV0A, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFV0A"), centFV0A, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FV0A, "Process with Run 3 FV0A estimator", false);

  void processRun3_FT0M(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::CentFT0Ms>::iterator const& col, aod::BCs const&)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator ft0m = initEstimator(col, "FT0");
    const float centFT0M = ft0m.getCentrality(col.multFT0M(), col.centFT0M());

    LOGF(debug, "centFT0M=%.0f", centFT0M);
    histos.fill(HIST("hCentFT0M"), centFT0M);
    histos.fill(HIST("hCentProfileFT0M"), centFT0M, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0M"), centFT0M, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0M, "Process with Run 3 FT0M estimator", false);

  void processRun3_FT0MOuterA(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::FITExtraMults, aod::CentFT0MOuterAs>::iterator const& col, aod::BCs const&)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator ft0mOuterA = initEstimator(col, "FT0MOuterA");
    const float centFT0MOuterA = ft0mOuterA.getCentrality(col.multFT0AOuter() + col.multFT0C(), col.centFT0MOuterA());

    LOGF(debug, "centFT0MOuterA=%.0f", centFT0MOuterA);
    histos.fill(HIST("hCentFT0MOuterA"), centFT0MOuterA);
    histos.fill(HIST("hCentProfileFT0MOuterA"), centFT0MOuterA, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0MOuterA"), centFT0MOuterA, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0MOuterA, "Process with Run 3 FT0M estimator", false);

  void processRun3_FT0A(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0As>::iterator const& col)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator ft0a = initEstimator(col, "FT0A");
    const float centFT0A = ft0a.getCentrality(col.multFT0A(), col.centFT0A());

    LOGF(debug, "centFT0A=%.0f", centFT0A);
    histos.fill(HIST("hCentFT0A"), centFT0A);
    histos.fill(HIST("hCentProfileFT0A"), centFT0A, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0A"), centFT0A, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0A, "Process with Run 3 FT0A estimator", false);

  void processRun3_FT0C(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::CentFT0Cs>::iterator const& col, aod::BCs const&)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator ft0c = initEstimator(col, "FT0C");
    const float centFT0C = ft0c.getCentrality(col.multFT0C(), col.centFT0C());

    LOGF(debug, "centFT0C=%.0f", centFT0C);
    histos.fill(HIST("hCentFT0C"), centFT0C);
    histos.fill(HIST("hCentProfileFT0C"), centFT0C, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0C"), centFT0C, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0C, "Process with Run 3 FT0C estimator", false);

  void processRun3_FT0CVar1(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::CentFT0CVariant1s>::iterator const& col, aod::BCs const&)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator ft0cVar1 = initEstimator(col, "FT0CVariant1");
    const float centFT0Cvar1 = ft0cVar1.getCentrality(col.multFT0C(), col.centFT0CVariant1());

    LOGF(debug, "centFT0Cvar1=%.0f", centFT0Cvar1);
    histos.fill(HIST("hCentFT0CVar1"), centFT0Cvar1);
    histos.fill(HIST("hCentProfileFT0CVar1"), centFT0Cvar1, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0CVar1"), centFT0Cvar1, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0CVar1, "Process with Run 3 FT0CVar1 estimator", false);

  void processRun3_FT0CVar2(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::CentFT0CVariant2s>::iterator const& col, aod::BCs const&)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator ft0cVar2 = initEstimator(col, "FT0CVariant2");
    const float centFT0Cvar2 = ft0cVar2.getCentrality(col.multFT0C(), col.centFT0CVariant2());

    LOGF(debug, "centFT0Cvar2=%.0f", centFT0Cvar2);
    histos.fill(HIST("hCentFT0CVar2"), centFT0Cvar2);
    histos.fill(HIST("hCentProfileFT0CVar2"), centFT0Cvar2, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0CVar2"), centFT0Cvar2, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0CVar2, "Process with Run 3 FT0CVar2 estimator", false);

  void processRun3_FDDM(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::CentFDDMs>::iterator const& col, aod::BCs const&)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator fddm = initEstimator(col, "FDDM");
    const float centFDDM = fddm.getCentrality(col.multFDDM(), col.centFDDM());

    LOGF(debug, "centFDDM=%.0f", centFDDM);
    histos.fill(HIST("hCentFDDM"), centFDDM);
    histos.fill(HIST("hCentProfileFDDM"), centFDDM, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFDDM"), centFDDM, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FDDM, "Process with Run 3 FDDM estimator", false);

  void processRun3_NTPV(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::CentNTPVs>::iterator const& col, aod::BCs const&)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator ntpv = initEstimator(col, "NTPV");
    const float centNTPV = ntpv.getCentrality(col.multNTracksPV(), col.centNTPV());

    LOGF(debug, "centNTPV=%.0f", centNTPV);
    histos.fill(HIST("hCentNTPV"), centNTPV);
    histos.fill(HIST("hCentProfileNTPV"), centNTPV, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentNTPV"), centNTPV, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_NTPV, "Process with Run 3 NTPV estimator", false);

  void processRun3_NGlobal(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::MultsGlobal, aod::CentNGlobals>::iterator const& col, aod::BCs const&)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator nGlo = initEstimator(col, "nGlo");
    const float centNGlo = nGlo.getCentrality(col.multNTracksGlobal(), col.centNGlobal());

    LOGF(debug, "centNGlo=%.0f", centNGlo);
    histos.fill(HIST("hCentNGlobal"), centNGlo);
    histos.fill(HIST("hCentProfileNGlobal"), centNGlo, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentNGlobal"), centNGlo, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_NGlobal, "Process with Run 3 NGlobal estimator", false);

  void processRun3_MFT(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::MFTMults, aod::CentMFTs>::iterator const& col, aod::BCs const&)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator mft = initEstimator(col, "MFT");
    const float centMFT = mft.getCentrality(col.mftNtracks(), col.centMFT());

    LOGF(debug, "centMFT=%.0f", centMFT);
    histos.fill(HIST("hCentMFT"), centMFT);
    histos.fill(HIST("hCentProfileMFT"), centMFT, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentMFT"), centMFT, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_MFT, "Process with Run 3 MFT estimator", false);

  void processRun3_FT0MAnchorCol(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::CentFT0MAnchorCols>::iterator const& col)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator ft0mAnchorCol = initEstimator(col, "FT0MAnchorCol");
    const float centFT0MAnchorCol = ft0mAnchorCol.getCentrality(col.multFT0M(), col.centFT0MAnchorCol());

    LOGF(debug, "centFT0MAnchorCol=%.0f", centFT0MAnchorCol);
    histos.fill(HIST("hCentFT0MAnchorCols"), centFT0MAnchorCol);
    histos.fill(HIST("hCentProfileFT0MAnchorCols"), centFT0MAnchorCol, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0MAnchorCols"), centFT0MAnchorCol, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0MAnchorCol, "Process with Run 3 FT0MAnchorCol estimator", false);

  void processRun3_FT0MAnchorBC(soa::Join<aod::Collisions, aod::EvSels, aod::MultsRun3, aod::CentFT0MAnchorBCs>::iterator const& col)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    Estimator ft0mAnchorBc = initEstimator(col, "FT0MAnchorBc");
    const float centFT0MAnchorBc = ft0mAnchorBc.getCentrality(col.multFT0M(), col.centFT0MAnchorBC());

    LOGF(debug, "centFT0MAnchorBc=%.0f", centFT0MAnchorBc);
    histos.fill(HIST("hCentFT0MAnchorBCs"), centFT0MAnchorBc);
    histos.fill(HIST("hCentProfileFT0MAnchorBCs"), centFT0MAnchorBc, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0MAnchorBCs"), centFT0MAnchorBc, col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processRun3_FT0MAnchorBC, "Process with Run 3 FT0MAnchorBC estimator", false);

  void processMonteCarloRun3_FV0A(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFV0As> const& collisions,
                                  soa::Join<aod::McCollisions, aod::MultMCExtras> const& mcCollisions,
                                  soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    for (auto const& col : collisions) {
      if (!isCollisionAccepted(col)) {
        continue;
      }

      Estimator fv0a = initEstimator(col, "FV0");
      const float centFV0A = fv0a.getCentrality(col.multFV0A(), col.centFV0A());

      LOGF(debug, "centFV0A=%.0f", centFV0A);
      histos.fill(HIST("hCentFV0A"), centFV0A);
      histos.fill(HIST("hCentProfileFV0A"), centFV0A, col.multNTracksPVetaHalf());
      histos.fill(HIST("hMultEta05VsCentFV0A"), centFV0A, col.multNTracksPVetaHalf());
      if (!loopOverMcCollisionsForMcHist) {
        const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
        histos.fill(HIST("hMultEta05VsGenMultFV0A"), mcCol.multMCFV0A(), col.multNTracksPVetaHalf());
        histos.fill(HIST("hGenMultEta05VsCentralityFV0A"), centFV0A, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFV0A"), centFV0A, mcCol.multMCFV0A());
      }
    }

    if (loopOverMcCollisionsForMcHist) {
      for (auto const& mcCol : mcCollisions) {
        auto groupedCollisions = collisions.sliceBy(perMcCollisionFV0A, mcCol.globalIndex());

        // Check if there is at least one of the reconstructed collisions associated to this MC collision
        // If so, we consider it
        int biggestNContribs = -1;
        int nContribsInEta05 = -1;
        float centrality = 100.5f;
        for (auto const& col : groupedCollisions) {
          if (!isCollisionAccepted(col)) {
            continue;
          }

          Estimator fv0a = initEstimator(col, "FV0");
          const float centFV0A = fv0a.getCentrality(col.multFV0A(), col.centFV0A());

          if (biggestNContribs < col.multPVTotalContributors()) {
            biggestNContribs = col.multPVTotalContributors();
            nContribsInEta05 = col.multNTracksPVetaHalf();
            centrality = centFV0A;
          }
        }

        histos.fill(HIST("hMultEta05VsGenMultFV0A"), mcCol.multMCFV0A(), nContribsInEta05);
        histos.fill(HIST("hGenMultEta05VsCentralityFV0A"), centrality, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFV0A"), centrality, mcCol.multMCFV0A());
      }
    }
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FV0A, "Process with Run 3 FV0A estimator", false);

  void processMonteCarloRun3_FT0M(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFT0Ms> const& collisions,
                                  soa::Join<aod::McCollisions, aod::MultMCExtras> const& mcCollisions,
                                  soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    for (auto const& col : collisions) {
      if (!isCollisionAccepted(col)) {
        return;
      }

      Estimator ft0m = initEstimator(col, "FT0");
      const float centFT0M = ft0m.getCentrality(col.multFT0M(), col.centFT0M());

      LOGF(debug, "centFT0M=%.0f", centFT0M);
      histos.fill(HIST("hCentFT0M"), centFT0M);
      histos.fill(HIST("hCentProfileFT0M"), centFT0M, col.multNTracksPVetaHalf());
      histos.fill(HIST("hMultEta05VsCentFT0M"), centFT0M, col.multNTracksPVetaHalf());
      if (!loopOverMcCollisionsForMcHist) {
        const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
        histos.fill(HIST("hMultEta05VsGenMultFT0M"), mcCol.multMCFT0A() + mcCol.multMCFT0C(), col.multNTracksPVetaHalf());
        histos.fill(HIST("hGenMultEta05VsCentralityFT0M"), centFT0M, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFT0M"), centFT0M, mcCol.multMCFT0A() + mcCol.multMCFT0C());
      }
    }

    if (loopOverMcCollisionsForMcHist) {
      for (auto const& mcCol : mcCollisions) {
        auto groupedCollisions = collisions.sliceBy(perMcCollisionFT0M, mcCol.globalIndex());

        // Check if there is at least one of the reconstructed collisions associated to this MC collision
        // If so, we consider it
        int biggestNContribs = -1;
        int nContribsInEta05 = -1;
        float centrality = 100.5f;
        for (auto const& col : groupedCollisions) {
          if (!isCollisionAccepted(col)) {
            continue;
          }

          Estimator ft0m = initEstimator(col, "FT0");
          const float centFT0M = ft0m.getCentrality(col.multFT0M(), col.centFT0M());

          if (biggestNContribs < col.multPVTotalContributors()) {
            biggestNContribs = col.multPVTotalContributors();
            nContribsInEta05 = col.multNTracksPVetaHalf();
            centrality = centFT0M;
          }
        }

        histos.fill(HIST("hMultEta05VsGenMultFT0M"), mcCol.multMCFT0A() + mcCol.multMCFT0C(), nContribsInEta05);
        histos.fill(HIST("hGenMultEta05VsCentralityFT0M"), centrality, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFT0M"), centrality, mcCol.multMCFT0A() + mcCol.multMCFT0C());
      }
    }
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FT0M, "Process with Run 3 FT0M estimator", false);

  void processMonteCarloRun3_FT0A(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFT0As> const& collisions,
                                  soa::Join<aod::McCollisions, aod::MultMCExtras> const& mcCollisions,
                                  soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    for (auto const& col : collisions) {
      if (!isCollisionAccepted(col)) {
        return;
      }

      Estimator ft0a = initEstimator(col, "FT0A");
      const float centFT0A = ft0a.getCentrality(col.multFT0A(), col.centFT0A());

      LOGF(debug, "centFT0M=%.0f", centFT0A);
      histos.fill(HIST("hCentFT0A"), centFT0A);
      histos.fill(HIST("hCentProfileFT0A"), centFT0A, col.multNTracksPVetaHalf());
      histos.fill(HIST("hMultEta05VsCentFT0A"), centFT0A, col.multNTracksPVetaHalf());
      if (!loopOverMcCollisionsForMcHist) {
        const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
        histos.fill(HIST("hMultEta05VsGenMultFT0A"), mcCol.multMCFT0A(), col.multNTracksPVetaHalf());
        histos.fill(HIST("hGenMultEta05VsCentralityFT0A"), centFT0A, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFT0A"), centFT0A, mcCol.multMCFT0A());
      }
    }

    if (loopOverMcCollisionsForMcHist) {
      for (auto const& mcCol : mcCollisions) {
        auto groupedCollisions = collisions.sliceBy(perMcCollisionFT0A, mcCol.globalIndex());

        // Check if there is at least one of the reconstructed collisions associated to this MC collision
        // If so, we consider it
        int biggestNContribs = -1;
        int nContribsInEta05 = -1;
        float centrality = 100.5f;
        for (auto const& col : groupedCollisions) {
          if (!isCollisionAccepted(col)) {
            continue;
          }

          Estimator ft0a = initEstimator(col, "FT0A");
          const float centFT0A = ft0a.getCentrality(col.multFT0A(), col.centFT0A());

          if (biggestNContribs < col.multPVTotalContributors()) {
            biggestNContribs = col.multPVTotalContributors();
            nContribsInEta05 = col.multNTracksPVetaHalf();
            centrality = centFT0A;
          }
        }

        histos.fill(HIST("hMultEta05VsGenMultFT0A"), mcCol.multMCFT0A(), nContribsInEta05);
        histos.fill(HIST("hGenMultEta05VsCentralityFT0A"), centrality, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFT0A"), centrality, mcCol.multMCFT0A());
      }
    }
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FT0A, "Process with Run 3 FT0A estimator", false);

  void processMonteCarloRun3_FT0C(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFT0Cs> const& collisions,
                                  soa::Join<aod::McCollisions, aod::MultMCExtras> const& mcCollisions,
                                  soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    for (auto const& col : collisions) {
      if (!isCollisionAccepted(col)) {
        return;
      }

      Estimator ft0c = initEstimator(col, "FT0C");
      const float centFT0C = ft0c.getCentrality(col.multFT0C(), col.centFT0C());

      LOGF(debug, "centFT0C=%.0f", centFT0C);
      histos.fill(HIST("hCentFT0C"), centFT0C);
      histos.fill(HIST("hCentProfileFT0C"), centFT0C, col.multNTracksPVetaHalf());
      histos.fill(HIST("hMultEta05VsCentFT0C"), centFT0C, col.multNTracksPVetaHalf());
      if (!loopOverMcCollisionsForMcHist) {
        const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
        histos.fill(HIST("hMultEta05VsGenMultFT0C"), mcCol.multMCFT0C(), col.multNTracksPVetaHalf());
        histos.fill(HIST("hGenMultEta05VsCentralityFT0C"), centFT0C, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFT0C"), centFT0C, mcCol.multMCFT0C());
      }
    }

    if (loopOverMcCollisionsForMcHist) {
      for (auto const& mcCol : mcCollisions) {
        auto groupedCollisions = collisions.sliceBy(perMcCollisionFT0C, mcCol.globalIndex());

        // Check if there is at least one of the reconstructed collisions associated to this MC collision
        // If so, we consider it
        int biggestNContribs = -1;
        int nContribsInEta05 = -1;
        float centrality = 100.5f;
        for (auto const& col : groupedCollisions) {
          if (!isCollisionAccepted(col)) {
            continue;
          }

          Estimator ft0c = initEstimator(col, "FT0C");
          const float centFT0C = ft0c.getCentrality(col.multFT0C(), col.centFT0C());

          if (biggestNContribs < col.multPVTotalContributors()) {
            biggestNContribs = col.multPVTotalContributors();
            nContribsInEta05 = col.multNTracksPVetaHalf();
            centrality = centFT0C;
          }
        }

        histos.fill(HIST("hMultEta05VsGenMultFT0C"), mcCol.multMCFT0C(), nContribsInEta05);
        histos.fill(HIST("hGenMultEta05VsCentralityFT0C"), centrality, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFT0C"), centrality, mcCol.multMCFT0C());
      }
    }
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FT0C, "Process with Run 3 FT0C estimator", false);

  void processMonteCarloRun3_FT0CVar1(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFT0CVariant1s> const& collisions,
                                      soa::Join<aod::McCollisions, aod::MultMCExtras> const& mcCollisions,
                                      soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    for (auto const& col : collisions) {
      if (!isCollisionAccepted(col)) {
        return;
      }

      Estimator ft0cVar1 = initEstimator(col, "FT0CVariant1");
      const float centFT0Cvar1 = ft0cVar1.getCentrality(col.multFT0C(), col.centFT0CVariant1());

      LOGF(debug, "centFT0Cvar1=%.0f", centFT0Cvar1);
      histos.fill(HIST("hCentFT0CVar1"), centFT0Cvar1);
      histos.fill(HIST("hCentProfileFT0CVar1"), centFT0Cvar1, col.multNTracksPVetaHalf());
      histos.fill(HIST("hMultEta05VsCentFT0CVar1"), centFT0Cvar1, col.multNTracksPVetaHalf());
      if (!loopOverMcCollisionsForMcHist) {
        const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
        histos.fill(HIST("hMultEta05VsGenMultFT0CVar1"), mcCol.multMCFT0C(), col.multNTracksPVetaHalf());
        histos.fill(HIST("hGenMultEta05VsCentralityFT0CVar1"), centFT0Cvar1, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFT0CVar1"), centFT0Cvar1, mcCol.multMCFT0C());
      }
    }

    if (loopOverMcCollisionsForMcHist) {
      for (auto const& mcCol : mcCollisions) {
        auto groupedCollisions = collisions.sliceBy(perMcCollisionFT0CVar1, mcCol.globalIndex());

        // Check if there is at least one of the reconstructed collisions associated to this MC collision
        // If so, we consider it
        int biggestNContribs = -1;
        int nContribsInEta05 = -1;
        float centrality = 100.5f;
        for (auto const& col : groupedCollisions) {
          if (!isCollisionAccepted(col)) {
            continue;
          }

          Estimator ft0cVar1 = initEstimator(col, "FT0CVariant1");
          const float centFT0Cvar1 = ft0cVar1.getCentrality(col.multFT0C(), col.centFT0CVariant1());

          if (biggestNContribs < col.multPVTotalContributors()) {
            biggestNContribs = col.multPVTotalContributors();
            nContribsInEta05 = col.multNTracksPVetaHalf();
            centrality = centFT0Cvar1;
          }
        }

        histos.fill(HIST("hMultEta05VsGenMultFT0CVar1"), mcCol.multMCFT0C(), nContribsInEta05);
        histos.fill(HIST("hGenMultEta05VsCentralityFT0CVar1"), centrality, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFT0CVar1"), centrality, mcCol.multMCFT0C());
      }
    }
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FT0CVar1, "Process with Run 3 FT0CVar1 estimator", false);

  void processMonteCarloRun3_FT0CVar2(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFT0CVariant2s> const& collisions,
                                      soa::Join<aod::McCollisions, aod::MultMCExtras> const& mcCollisions,
                                      soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    for (auto const& col : collisions) {
      if (!isCollisionAccepted(col)) {
        return;
      }

      Estimator ft0cVar2 = initEstimator(col, "FT0CVariant2");
      const float centFT0Cvar2 = ft0cVar2.getCentrality(col.multFT0C(), col.centFT0CVariant2());

      LOGF(debug, "centFT0Cvar2=%.0f", centFT0Cvar2);
      histos.fill(HIST("hCentFT0CVar2"), centFT0Cvar2);
      histos.fill(HIST("hCentProfileFT0CVar2"), centFT0Cvar2, col.multNTracksPVetaHalf());
      histos.fill(HIST("hMultEta05VsCentFT0CVar2"), centFT0Cvar2, col.multNTracksPVetaHalf());
      if (!loopOverMcCollisionsForMcHist) {
        const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
        histos.fill(HIST("hMultEta05VsGenMultFT0CVar2"), mcCol.multMCFT0C(), col.multNTracksPVetaHalf());
        histos.fill(HIST("hGenMultEta05VsCentralityFT0CVar2"), centFT0Cvar2, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFT0CVar2"), centFT0Cvar2, mcCol.multMCFT0C());
      }
    }

    if (loopOverMcCollisionsForMcHist) {
      for (auto const& mcCol : mcCollisions) {
        auto groupedCollisions = collisions.sliceBy(perMcCollisionFT0CVar2, mcCol.globalIndex());

        // Check if there is at least one of the reconstructed collisions associated to this MC collision
        // If so, we consider it
        int biggestNContribs = -1;
        int nContribsInEta05 = -1;
        float centrality = 100.5f;
        for (auto const& col : groupedCollisions) {
          if (!isCollisionAccepted(col)) {
            continue;
          }

          Estimator ft0cVar2 = initEstimator(col, "FT0CVariant2");
          const float centFT0Cvar2 = ft0cVar2.getCentrality(col.multFT0C(), col.centFT0CVariant2());

          if (biggestNContribs < col.multPVTotalContributors()) {
            biggestNContribs = col.multPVTotalContributors();
            nContribsInEta05 = col.multNTracksPVetaHalf();
            centrality = centFT0Cvar2;
          }
        }

        histos.fill(HIST("hMultEta05VsGenMultFT0CVar2"), mcCol.multMCFT0C(), nContribsInEta05);
        histos.fill(HIST("hGenMultEta05VsCentralityFT0CVar2"), centrality, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFT0CVar2"), centrality, mcCol.multMCFT0C());
      }
    }
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FT0CVar2, "Process with Run 3 FT0CVar2 estimator", false);

  void processMonteCarloRun3_FDDM(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentFDDMs> const& collisions,
                                  soa::Join<aod::McCollisions, aod::MultMCExtras> const& mcCollisions,
                                  soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    for (auto const& col : collisions) {
      if (!isCollisionAccepted(col)) {
        return;
      }

      Estimator fddm = initEstimator(col, "FDDM");
      const float centFDDM = fddm.getCentrality(col.multFDDM(), col.centFDDM());

      LOGF(debug, "centFDDM=%.0f", centFDDM);
      histos.fill(HIST("hCentFDDM"), centFDDM);
      histos.fill(HIST("hCentProfileFDDM"), centFDDM, col.multNTracksPVetaHalf());
      histos.fill(HIST("hMultEta05VsCentFDDM"), centFDDM, col.multNTracksPVetaHalf());
      if (!loopOverMcCollisionsForMcHist) {
        const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
        histos.fill(HIST("hMultEta05VsGenMultFDDM"), mcCol.multMCFDDA() + mcCol.multMCFDDC(), col.multNTracksPVetaHalf());
        histos.fill(HIST("hGenMultEta05VsCentralityFDDM"), centFDDM, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFDDM"), centFDDM, mcCol.multMCFDDA() + mcCol.multMCFDDC());
      }
    }

    if (loopOverMcCollisionsForMcHist) {
      for (auto const& mcCol : mcCollisions) {
        auto groupedCollisions = collisions.sliceBy(perMcCollisionFDDM, mcCol.globalIndex());

        // Check if there is at least one of the reconstructed collisions associated to this MC collision
        // If so, we consider it
        int biggestNContribs = -1;
        int nContribsInEta05 = -1;
        float centrality = 100.5f;
        for (auto const& col : groupedCollisions) {
          if (!isCollisionAccepted(col)) {
            continue;
          }

          Estimator fddm = initEstimator(col, "FDDM");
          const float centFDDM = fddm.getCentrality(col.multFDDM(), col.centFDDM());

          if (biggestNContribs < col.multPVTotalContributors()) {
            biggestNContribs = col.multPVTotalContributors();
            nContribsInEta05 = col.multNTracksPVetaHalf();
            centrality = centFDDM;
          }
        }

        histos.fill(HIST("hMultEta05VsGenMultFDDM"), mcCol.multMCFDDA() + mcCol.multMCFDDC(), nContribsInEta05);
        histos.fill(HIST("hGenMultEta05VsCentralityFDDM"), centrality, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityFDDM"), centrality, mcCol.multMCFDDA() + mcCol.multMCFDDC());
      }
    }
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FDDM, "Process with Run 3 FDDM estimator", false);

  void processMonteCarloRun3_NTPV(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsExtra, aod::CentNTPVs> const& collisions,
                                  soa::Join<aod::McCollisions, aod::MultMCExtras> const& mcCollisions,
                                  soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    for (auto const& col : collisions) {
      if (!isCollisionAccepted(col)) {
        return;
      }

      Estimator ntpv = initEstimator(col, "NTPV");
      const float centNTPV = ntpv.getCentrality(col.multNTracksPV(), col.centNTPV());

      histos.fill(HIST("hCentNTPV"), centNTPV);
      histos.fill(HIST("hCentProfileNTPV"), centNTPV, col.multNTracksPVetaHalf());
      histos.fill(HIST("hMultEta05VsCentNTPV"), centNTPV, col.multNTracksPVetaHalf());
      if (!loopOverMcCollisionsForMcHist) {
        const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
        histos.fill(HIST("hMultEta05VsGenMultNTPV"), mcCol.multMCNParticlesEta08(), col.multNTracksPVetaHalf());
        histos.fill(HIST("hGenMultEta05VsCentralityNTPV"), centNTPV, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityNTPV"), centNTPV, mcCol.multMCNParticlesEta08());
      }
    }

    if (loopOverMcCollisionsForMcHist) {
      for (auto const& mcCol : mcCollisions) {
        auto groupedCollisions = collisions.sliceBy(perMcCollisionNTPV, mcCol.globalIndex());

        // Check if there is at least one of the reconstructed collisions associated to this MC collision
        // If so, we consider it
        int biggestNContribs = -1;
        int nContribsInEta05 = -1;
        float centrality = 100.5f;
        for (auto const& col : groupedCollisions) {
          if (!isCollisionAccepted(col)) {
            continue;
          }

          Estimator ntpv = initEstimator(col, "NTPV");
          const float centNTPV = ntpv.getCentrality(col.multNTracksPV(), col.centNTPV());

          if (biggestNContribs < col.multPVTotalContributors()) {
            biggestNContribs = col.multPVTotalContributors();
            nContribsInEta05 = col.multNTracksPVetaHalf();
            centrality = centNTPV;
          }
        }

        histos.fill(HIST("hMultEta05VsGenMultNTPV"), mcCol.multMCNParticlesEta08(), nContribsInEta05);
        histos.fill(HIST("hGenMultEta05VsCentralityNTPV"), centrality, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityNTPV"), centrality, mcCol.multMCNParticlesEta08());
      }
    }
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_NTPV, "Process with Run 3 NTPV estimator", false);

  void processMonteCarloRun3_NGlobal(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MultsGlobal, aod::MultsExtra, aod::CentNGlobals> const& collisions,
                                     soa::Join<aod::McCollisions, aod::MultMCExtras> const& mcCollisions,
                                     soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    for (auto const& col : collisions) {
      if (!isCollisionAccepted(col)) {
        return;
      }

      Estimator nGlo = initEstimator(col, "nGlo");
      const float centNGlo = nGlo.getCentrality(col.multNTracksGlobal(), col.centNGlobal());

      LOGF(debug, "centNGlo=%.0f", centNGlo);
      histos.fill(HIST("hCentNGlobal"), centNGlo);
      histos.fill(HIST("hCentProfileNGlobal"), centNGlo, col.multNTracksPVetaHalf());
      histos.fill(HIST("hMultEta05VsCentNGlobal"), centNGlo, col.multNTracksPVetaHalf());
      if (!loopOverMcCollisionsForMcHist) {
        const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
        histos.fill(HIST("hMultEta05VsGenMultNGlobal"), mcCol.multMCNParticlesEta08(), col.multNTracksPVetaHalf());
        histos.fill(HIST("hGenMultEta05VsCentralityNGlobal"), centNGlo, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityNGlobal"), centNGlo, mcCol.multMCNParticlesEta08());
      }
    }

    if (loopOverMcCollisionsForMcHist) {
      for (auto const& mcCol : mcCollisions) {
        auto groupedCollisions = collisions.sliceBy(perMcCollisionNGlobal, mcCol.globalIndex());

        // Check if there is at least one of the reconstructed collisions associated to this MC collision
        // If so, we consider it
        int biggestNContribs = -1;
        int nContribsInEta05 = -1;
        float centrality = 100.5f;
        for (auto const& col : groupedCollisions) {
          if (!isCollisionAccepted(col)) {
            continue;
          }

          Estimator nGlo = initEstimator(col, "nGlo");
          const float centNGlo = nGlo.getCentrality(col.multNTracksGlobal(), col.centNGlobal());

          if (biggestNContribs < col.multPVTotalContributors()) {
            biggestNContribs = col.multPVTotalContributors();
            nContribsInEta05 = col.multNTracksPVetaHalf();
            centrality = centNGlo;
          }
        }

        histos.fill(HIST("hMultEta05VsGenMultNGlobal"), mcCol.multMCNParticlesEta08(), nContribsInEta05);
        histos.fill(HIST("hGenMultEta05VsCentralityNGlobal"), centrality, mcCol.multMCNParticlesEta05());
        histos.fill(HIST("hGenMultVsCentralityNGlobal"), centrality, mcCol.multMCNParticlesEta08());
      }
    }
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_NGlobal, "Process with Run 3 NGlobal estimator", false);

  void processMonteCarloRun3_MFT(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::MFTMults, aod::MultsExtra, aod::CentMFTs> const& collisions,
                                 soa::Join<aod::McCollisions, aod::MultMCExtras> const& mcCollisions,
                                 soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    for (auto const& col : collisions) {
      if (!isCollisionAccepted(col)) {
        return;
      }

      Estimator mft = initEstimator(col, "MFT");
      const float centMFT = mft.getCentrality(col.mftNtracks(), col.centMFT());

      LOGF(debug, "centMFT=%.0f", centMFT);
      histos.fill(HIST("hCentMFT"), centMFT);
      histos.fill(HIST("hCentProfileMFT"), centMFT, col.multNTracksPVetaHalf());
      histos.fill(HIST("hMultEta05VsCentMFT"), centMFT, col.multNTracksPVetaHalf());
      if (!loopOverMcCollisionsForMcHist) {
        const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
        // histos.fill(HIST("hMultEta05VsGenMultMFT"), mcCol.multMCMFT(), col.multNTracksPVetaHalf()); // FIXME: uncomment when MC MFT mult is added in aod::MultMCExtras
        histos.fill(HIST("hGenMultEta05VsCentralityMFT"), centMFT, mcCol.multMCNParticlesEta05());
        // histos.fill(HIST("hGenMultVsCentralityMFT"), centMFT, mcCol.multMCMFT()); // FIXME: uncomment when MC MFT mult is added in aod::MultMCExtras
      }
    }

    if (loopOverMcCollisionsForMcHist) {
      for (auto const& mcCol : mcCollisions) {
        auto groupedCollisions = collisions.sliceBy(perMcCollisionMFT, mcCol.globalIndex());

        // Check if there is at least one of the reconstructed collisions associated to this MC collision
        // If so, we consider it
        int biggestNContribs = -1;
        // int nContribsInEta05 = -1;
        float centrality = 100.5f;
        for (auto const& col : groupedCollisions) {
          if (!isCollisionAccepted(col)) {
            continue;
          }

          Estimator mft = initEstimator(col, "MFT");
          const float centMFT = mft.getCentrality(col.mftNtracks(), col.centMFT());

          if (biggestNContribs < col.multPVTotalContributors()) {
            biggestNContribs = col.multPVTotalContributors();
            // nContribsInEta05 = col.multNTracksPVetaHalf();
            centrality = centMFT;
          }
        }

        // histos.fill(HIST("hMultEta05VsGenMultMFT"), mcCol.multMCMFT(), nContribsInEta05); // FIXME: uncomment when MC MFT mult is added in aod::MultMCExtras
        histos.fill(HIST("hGenMultEta05VsCentralityMFT"), centrality, mcCol.multMCNParticlesEta05());
        // histos.fill(HIST("hGenMultVsCentralityMFT"), col.centMFT(), mcCol.multMCMFT()); // FIXME: uncomment when MC MFT mult is added in aod::MultMCExtras
      }
    }
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_MFT, "Process with Run 3 MFT estimator", false);

  void processMonteCarloRun3_FT0MAnchorCol(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::CentFT0MAnchorCols>::iterator const& col,
                                           soa::Join<aod::McCollisions, aod::MultMCExtras> const& /*mcCollisions*/,
                                           soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }

    const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
    Estimator ft0mAnchorCol = initEstimator(col, "FT0MAnchorCol");
    const float centFT0MAnchorCol = ft0mAnchorCol.getCentrality(col.multFT0M(), col.centFT0MAnchorCol());

    LOGF(debug, "centFT0MAnchorCol=%.0f", centFT0MAnchorCol);
    histos.fill(HIST("hCentFT0MAnchorCols"), centFT0MAnchorCol);
    histos.fill(HIST("hCentProfileFT0MAnchorCols"), centFT0MAnchorCol, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0MAnchorCols"), centFT0MAnchorCol, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsGenMultFT0MAnchorCols"), mcCol.multMCFT0A() + mcCol.multMCFT0C(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FT0MAnchorCol, "Process with Run 3 FT0MAnchorCol estimator", false);

  void processMonteCarloRun3_FT0MAnchorBC(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::MultsRun3, aod::CentFT0MAnchorBCs>::iterator const& col,
                                          soa::Join<aod::McCollisions, aod::MultMCExtras> const& /*mcCollisions*/,
                                          soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/)
  {
    if (!isCollisionAccepted(col)) {
      return;
    }
    const auto& mcCol = col.mcCollision_as<soa::Join<aod::McCollisions, aod::MultMCExtras>>();
    Estimator ft0mAnchorBc = initEstimator(col, "FT0MAnchorBc");
    const float centFT0MAnchorBc = ft0mAnchorBc.getCentrality(col.multFT0M(), col.centFT0MAnchorBC());

    LOGF(debug, "centFT0MAnchorBc=%.0f", centFT0MAnchorBc);
    histos.fill(HIST("hCentFT0MAnchorBCs"), centFT0MAnchorBc);
    histos.fill(HIST("hCentProfileFT0MAnchorBCs"), centFT0MAnchorBc, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsCentFT0MAnchorBCs"), centFT0MAnchorBc, col.multNTracksPVetaHalf());
    histos.fill(HIST("hMultEta05VsGenMultFT0MAnchorBCs"), mcCol.multMCFT0A() + mcCol.multMCFT0C(), col.multNTracksPVetaHalf());
  }
  PROCESS_SWITCH(CentralityQa, processMonteCarloRun3_FT0MAnchorBC, "Process with Run 3 FT0MAnchorBC estimator", false);

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
  PROCESS_SWITCH(CentralityQa, processBunchCrossings, "Process with Run 3 BC table", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CentralityQa>(cfgc)};
}
