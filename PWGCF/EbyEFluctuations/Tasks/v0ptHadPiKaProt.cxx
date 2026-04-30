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

/// \file v0ptHadPiKaProt.cxx
/// \brief Task for analyzing v0(pT) of inclusive hadrons, pions, kaons, and, protons
/// \author Swati Saha

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include "TDatabasePDG.h"
#include <TComplex.h>
#include <TF1.h>
#include <TH2.h>
#include <THn.h>
#include <TList.h>
#include <TMath.h>
#include <TPDGCode.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom3.h>
#include <TString.h>

#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static constexpr float LongArrayFloat[3][20] = {{1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {2.1, 2.2, 2.3, -2.1, -2.2, -2.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {3.1, 3.2, 3.3, -3.1, -3.2, -3.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}};

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct V0ptHadPiKaProt {

  // ITS response
  o2::aod::ITSResponse itsResponse;
  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "https://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Users/s/swati/PhiWeight", "CCDB path to ccdb object containing phi weight in a 3D histogram"};
  Configurable<int64_t> ccdbNoLaterThanPtEff{"ccdbNoLaterThanPtEff", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbPathPtEff{"ccdbPathPtEff", "Users/s/swati/EfficiencyWeight", "CCDB path to ccdb object containing pt-dependent efficiency"};

  enum Particles {
    PIONS = 0,
    KAONS,
    PROTONS
  };
  enum ParticleNsigma {
    kPionUpCut = 0,
    kKaonUpCut,
    kProtonUpCut,
    kPionLowCut,
    kKaonLowCut,
    kProtonLowCut
  };
  enum DetectorType {
    kTPC = 0,
    kTOF,
    kITS
  };
  enum CentralityEstimator {
    kFT0C = 0,
    kFT0A,
    kFT0M,
    kFV0A
  };

  Configurable<bool> cfgIsMC{"cfgIsMC", false, "Run MC"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutTpcChi2NCl{"cfgCutTpcChi2NCl", 2.5f, "Maximum TPCchi2NCl"};
  Configurable<float> cfgCutItsChi2NCl{"cfgCutItsChi2NCl", 36.0f, "Maximum ITSchi2NCl"};
  Configurable<float> cfgCutTrackDcaZ{"cfgCutTrackDcaZ", 2.0f, "Maximum DcaZ"};
  Configurable<int> cfgITScluster{"cfgITScluster", 1, "Minimum Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 80, "Minimum Number of TPC cluster"};
  Configurable<int> cfgTPCnCrossedRows{"cfgTPCnCrossedRows", 70, "Minimum Number of TPC crossed-rows"};
  Configurable<float> cfgCutPtUpperTPC{"cfgCutPtUpperTPC", 0.6f, "Upper pT cut for PID using TPC only"};
  Configurable<float> cfgnSigmaOtherParticles{"cfgnSigmaOtherParticles", 3.0f, "PID nSigma cut to remove other particles (default:3)"};
  Configurable<float> cfgnSigmaCutTPC{"cfgnSigmaCutTPC", 2.0f, "PID nSigma cut for TPC"};
  Configurable<float> cfgnSigmaCutTOF{"cfgnSigmaCutTOF", 2.0f, "PID nSigma cut for TOF"};
  Configurable<bool> cfgUseNewSeperationPid{"cfgUseNewSeperationPid", true, "Use seperation based PID cuts (NEW)"};
  Configurable<float> cfgnSigmaCutTPCHigherPt{"cfgnSigmaCutTPCHigherPt", 2.0f, "PID nSigma cut for TPC at higher pt"};
  Configurable<float> cfgnSigmaCutTOFHigherPt{"cfgnSigmaCutTOFHigherPt", 2.0f, "PID nSigma cut for TOF at higher pt"};
  Configurable<float> cfgnSigmaSeperationCut{"cfgnSigmaSeperationCut", 3.5f, "PID nSigma of other species must be greater than the vale"};
  Configurable<float> cfgnSigmaCutCombTPCTOF{"cfgnSigmaCutCombTPCTOF", 2.0f, "PID nSigma combined cut for TPC and TOF"};
  ConfigurableAxis nchAxis{"nchAxis", {5000, 0.5, 5000.5}, ""};
  ConfigurableAxis centAxis{"centAxis", {90, 0., 90.}, "Centrality/Multiplicity percentile bining"};
  ConfigurableAxis nchAxis1{"nchAxis1", {500, 0.5, 500.5}, "Axis for multiplicity of GlobalTracks/PVTracks"};
  ConfigurableAxis nchAxis2{"nchAxis2", {1000, 0.5, 30000.5}, "Axis for multiplicity of FT0A/FT0C/FV0A"};
  ConfigurableAxis nchAxis3{"nchAxis3", {1000, 0.5, 100000.5}, "Axis for multiplicity of FT0A/FT0C/FV0A"};
  Configurable<float> cfgCutPtLower{"cfgCutPtLower", 0.2f, "Lower pT cut"};
  Configurable<float> cfgCutPtLowerProt{"cfgCutPtLowerProt", 0.2f, "Lower pT cut"};
  Configurable<float> cfgCutPtUpper{"cfgCutPtUpper", 10.0f, "Higher pT cut for inclusive hadron analysis"};
  Configurable<float> cfgCutPtUpperPID{"cfgCutPtUpperPID", 6.0f, "Higher pT cut for identified particle analysis"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "absolute Eta cut"};
  Configurable<float> cfgCutEtaLeft{"cfgCutEtaLeft", -0.4f, "Left end of eta gap"};
  Configurable<float> cfgCutEtaRight{"cfgCutEtaRight", 0.4f, "Right end of eta gap"};
  Configurable<int> cfgNSubsample{"cfgNSubsample", 20, "Number of subsamples"};
  Configurable<int> cfgCentralityChoice{"cfgCentralityChoice", 0, "Which centrality estimator? 0-->FT0C, 1-->FT0A, 2-->FT0M, 3-->FV0A"};
  Configurable<bool> cfgEvSelkNoSameBunchPileup{"cfgEvSelkNoSameBunchPileup", true, "Pileup removal"};
  Configurable<bool> cfgUseGoodITSLayerAllCut{"cfgUseGoodITSLayerAllCut", true, "Remove time interval with dead ITS zone"};
  Configurable<bool> cfgEvSelkNoITSROFrameBorder{"cfgEvSelkNoITSROFrameBorder", true, "ITSROFrame border event selection cut"};
  Configurable<bool> cfgEvSelkNoTimeFrameBorder{"cfgEvSelkNoTimeFrameBorder", true, "TimeFrame border event selection cut"};
  Configurable<bool> cfgEvSelUseGoodZvtxFT0vsPV{"cfgEvSelUseGoodZvtxFT0vsPV", true, "GoodZvertex and FT0 vs PV cut"};
  Configurable<bool> cfgEvSelUseOcuppancyTimeCut{"cfgEvSelUseOcuppancyTimeCut", true, "Occupancy Time pattern cut"};
  Configurable<bool> cfgEvSelSetOcuppancyRange{"cfgEvSelSetOcuppancyRange", true, "Use cut on occupancy range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "min. value of occupancy"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 3000, "max. value of occupancy"};

  Configurable<bool> cfgUseItsPID{"cfgUseItsPID", false, "Use ITS PID for particle identification"};
  Configurable<float> cfgPtCutTOF{"cfgPtCutTOF", 0.3f, "Minimum pt to use TOF N-sigma"};
  Configurable<LabeledArray<float>> nSigmas{"nSigmas", {LongArrayFloat[0], 3, 6, {"TPC", "TOF", "ITS"}, {"pos_pi", "pos_ka", "pos_pr", "neg_pi", "neg_ka", "neg_pr"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};
  Configurable<bool> cfgUseRun3V2PID{"cfgUseRun3V2PID", true, "True if PID cuts to be used are similar to Run3 v2 PID analysis"};
  Configurable<int> cfgNbinsV02pt{"cfgNbinsV02pt", 14, "No. of pT bins for v02(pT) analysis"};
  Configurable<float> cfgCutPtMaxForV02{"cfgCutPtMaxForV02", 3.0f, "Max. pT for v02(pT)"};
  Configurable<float> cfgCutEtaWindowB{"cfgCutEtaWindowB", 0.4f, "value of x in |eta|<x for window B"};
  Configurable<bool> cfgLoadPhiWeights{"cfgLoadPhiWeights", false, "Load phi weights from CCDB to take care of non-uniform acceptance"};
  Configurable<bool> cfgLoadPtEffWeights{"cfgLoadPtEffWeights", false, "Load pt-dependent efficiency weights from CCDB to take care of detector inefficiency"};
  Configurable<int> cfgMinNoOfParticles{"cfgMinNoOfParticles", 4, "Minimum no. of particles for calculating v02(pT)"};
  Configurable<int> cfgV02WeightedFill{"cfgV02WeightedFill", false, "Fill profiles related to v2 with multiplicity-based weights?"};

  // pT dep DCAxy and DCAz cuts
  Configurable<bool> cfgUsePtDepDCAxy{"cfgUsePtDepDCAxy", true, "Use pt-dependent DCAxy cut"};
  Configurable<bool> cfgUsePtDepDCAz{"cfgUsePtDepDCAz", true, "Use pt-dependent DCAz cut"};
  O2_DEFINE_CONFIGURABLE(cfgDCAxyFunc, std::string, "(0.0026+0.005/(x^1.01))", "Functional form of pt-dependent DCAxy cut");
  O2_DEFINE_CONFIGURABLE(cfgDCAzFunc, std::string, "(0.0026+0.005/(x^1.01))", "Functional form of pt-dependent DCAz cut");
  TF1* fPtDepDCAxy = nullptr;
  TF1* fPtDepDCAz = nullptr;

  O2_DEFINE_CONFIGURABLE(cfgUseSmallIonAdditionalEventCut, bool, true, "Use additional event cut on mult correlations for small ions")
  O2_DEFINE_CONFIGURABLE(cfgEvSelMultCorrelation, bool, true, "Multiplicity correlation cut")
  O2_DEFINE_CONFIGURABLE(cfgEvSelV0AT0ACut, bool, true, "V0A T0A 5 sigma cut")
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

  // Output
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosAnalysis{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::vector<std::vector<std::shared_ptr<TProfile2D>>> subSample;
  std::vector<std::vector<std::shared_ptr<TProfile2D>>> subSampleV02;
  std::vector<std::vector<std::shared_ptr<TProfile2D>>> subSampleV02_weighted;
  TRandom3* funRndm = new TRandom3(0);

  // Phi weight histograms initialization
  TH3D* hWeightPhiFunctionVzEtaPhi = nullptr;
  // Efficiency of diff. particle histograms initialization
  TH1D* hEffAllCharged = nullptr;
  TH1D* hEffPion = nullptr;
  TH1D* hEffKaon = nullptr;
  TH1D* hEffProton = nullptr;

  // Filter command***********
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtLower) && (aod::track::pt < cfgCutPtUpper) && (requireGlobalTrackInFilter()) && (aod::track::tpcChi2NCl < cfgCutTpcChi2NCl) && (aod::track::itsChi2NCl < cfgCutItsChi2NCl) && (nabs(aod::track::dcaZ) < cfgCutTrackDcaZ);

  // Filtering collisions and tracks***********
  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs, aod::Mults>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullEl, aod::pidTOFFullEl>>;

  using MyMCRecCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::Mults, aod::McCollisionLabels>>;
  using MyMCTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullEl, aod::pidTOFFullEl>>;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  std::array<float, 6> tofNsigmaCut;
  std::array<float, 6> itsNsigmaCut;
  std::array<float, 6> tpcNsigmaCut;

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
    // Get nSigma values of TPC, TOF, ITS for various particles from the matrix "nSigmas"
    tpcNsigmaCut[kPionUpCut] = nSigmas->getData()[kTPC][kPionUpCut];
    tpcNsigmaCut[kKaonUpCut] = nSigmas->getData()[kTPC][kKaonUpCut];
    tpcNsigmaCut[kProtonUpCut] = nSigmas->getData()[kTPC][kProtonUpCut];
    tpcNsigmaCut[kPionLowCut] = nSigmas->getData()[kTPC][kPionLowCut];
    tpcNsigmaCut[kKaonLowCut] = nSigmas->getData()[kTPC][kKaonLowCut];
    tpcNsigmaCut[kProtonLowCut] = nSigmas->getData()[kTPC][kProtonLowCut];

    tofNsigmaCut[kPionUpCut] = nSigmas->getData()[kTOF][kPionUpCut];
    tofNsigmaCut[kKaonUpCut] = nSigmas->getData()[kTOF][kKaonUpCut];
    tofNsigmaCut[kProtonUpCut] = nSigmas->getData()[kTOF][kProtonUpCut];
    tofNsigmaCut[kPionLowCut] = nSigmas->getData()[kTOF][kPionLowCut];
    tofNsigmaCut[kKaonLowCut] = nSigmas->getData()[kTOF][kKaonLowCut];
    tofNsigmaCut[kProtonLowCut] = nSigmas->getData()[kTOF][kProtonLowCut];

    itsNsigmaCut[kPionUpCut] = nSigmas->getData()[kITS][kPionUpCut];
    itsNsigmaCut[kKaonUpCut] = nSigmas->getData()[kITS][kKaonUpCut];
    itsNsigmaCut[kProtonUpCut] = nSigmas->getData()[kITS][kProtonUpCut];
    itsNsigmaCut[kPionLowCut] = nSigmas->getData()[kITS][kPionLowCut];
    itsNsigmaCut[kKaonLowCut] = nSigmas->getData()[kITS][kKaonLowCut];
    itsNsigmaCut[kProtonLowCut] = nSigmas->getData()[kITS][kProtonLowCut];

    // Loading phi weight histograms from CCDB
    if (cfgLoadPhiWeights) {

      // Accessing phi weight histograms
      ccdb->setURL(ccdbUrl.value);
      // Enabling object caching, otherwise each call goes to the CCDB server
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      // Not later than now, will be replaced by the value of the train creation
      // This avoids that users can replace objects **while** a train is running
      ccdb->setCreatedNotAfter(ccdbNoLaterThan.value);
      LOGF(info, "Getting object %s", ccdbPath.value.data());
      TList* lst = ccdb->getForTimeStamp<TList>(ccdbPath.value, ccdbNoLaterThan.value);
      hWeightPhiFunctionVzEtaPhi = reinterpret_cast<TH3D*>(lst->FindObject("hWeightPhiFunctionVzEtaPhi"));
      if (!hWeightPhiFunctionVzEtaPhi)
        LOGF(info, "FATAL!! could not get phi weights---------> check");
    }

    // Loading pT-dependent efficiency histograms from CCDB
    if (cfgLoadPtEffWeights) {

      ccdb->setURL(ccdbUrl.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(ccdbNoLaterThanPtEff.value);
      LOGF(info, "Getting object %s", ccdbPathPtEff.value.data());
      TList* lst = ccdb->getForTimeStamp<TList>(ccdbPathPtEff.value, ccdbNoLaterThanPtEff.value);
      hEffAllCharged = reinterpret_cast<TH1D*>(lst->FindObject("hEffAllCharged"));
      hEffPion = reinterpret_cast<TH1D*>(lst->FindObject("hEffPion"));
      hEffKaon = reinterpret_cast<TH1D*>(lst->FindObject("hEffKaon"));
      hEffProton = reinterpret_cast<TH1D*>(lst->FindObject("hEffProton"));
      if (!hEffAllCharged || !hEffPion || !hEffKaon || !hEffProton)
        LOGF(info, "FATAL!! could not get efficiency files---------> !!! check !!!");
    }

    // Define axes
    std::vector<double> ptBin = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
    AxisSpec ptAxis = {ptBin, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec noAxis = {1, 0, 1, "no axis"};
    AxisSpec nSigmaAxis = {200, -5.0, 5.0, "n#sigma"};

    // Add histograms to histogram manager (as in the output object of in AliPhysics)

    // QA hists
    histos.add("hEventStatData", "Data Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
    histos.add("hZvtx_after_sel", ";Z (cm)", kTH1F, {{240, -12, 12}});
    histos.add("hCentrality", ";centrality (%)", kTH1F, {{90, 0, 90}});
    // before selection
    histos.add("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_PVTracks_beforeSel", "", {HistType::kTH2D, {nchAxis1, nchAxis1}});
    histos.add("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_centFT0C_beforeSel", "", {HistType::kTH2D, {centAxis, nchAxis1}});
    histos.add("MultCorrelationPlots/BeforeSelection/His2D_PVTracks_centFT0C_beforeSel", "", {HistType::kTH2D, {centAxis, nchAxis1}});
    histos.add("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_V0ATracks_beforeSel", "", {HistType::kTH2D, {nchAxis3, nchAxis1}});
    histos.add("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_T0ATracks_beforeSel", "", {HistType::kTH2D, {nchAxis2, nchAxis1}});
    histos.add("MultCorrelationPlots/BeforeSelection/His2D_V0ATracks_T0CTracks_beforeSel", "", {HistType::kTH2D, {nchAxis2, nchAxis3}});
    // after selection
    if (cfgUseSmallIonAdditionalEventCut) {
      histos.add("MultCorrelationPlots/AfterSelection/His2D_globalTracks_PVTracks_afterSel", "", {HistType::kTH2D, {nchAxis1, nchAxis1}});
      histos.add("MultCorrelationPlots/AfterSelection/His2D_globalTracks_centFT0C_afterSel", "", {HistType::kTH2D, {centAxis, nchAxis1}});
      histos.add("MultCorrelationPlots/AfterSelection/His2D_PVTracks_centFT0C_afterSel", "", {HistType::kTH2D, {centAxis, nchAxis1}});
      histos.add("MultCorrelationPlots/AfterSelection/His2D_globalTracks_V0ATracks_afterSel", "", {HistType::kTH2D, {nchAxis3, nchAxis1}});
      histos.add("MultCorrelationPlots/AfterSelection/His2D_globalTracks_T0ATracks_afterSel", "", {HistType::kTH2D, {nchAxis2, nchAxis1}});
      histos.add("MultCorrelationPlots/AfterSelection/His2D_V0ATracks_T0CTracks_afterSel", "", {HistType::kTH2D, {nchAxis2, nchAxis3}});
    }

    histos.add("Hist2D_globalTracks_PVTracks", "", {HistType::kTH2D, {nchAxis, nchAxis}});
    histos.add("Hist2D_cent_nch", "", {HistType::kTH2D, {nchAxis, centAxis}});
    histos.add("hP", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.2, 4.}});
    histos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hPhi", ";#phi", kTH1F, {{100, 0., o2::constants::math::TwoPI}});
    histos.add("hEta", ";#eta", kTH1F, {{100, -2.01, 2.01}});
    histos.add("hDcaXY", ";#it{dca}_{XY}", kTH1F, {{1000, -0.5, 0.5}});
    histos.add("hDcaZ", ";#it{dca}_{Z}", kTH1F, {{1000, -0.5, 0.5}});
    histos.add("hMeanPt", "", kTProfile, {centAxis});

    // phi weight hist for non-uniform acceptance correction
    histos.add("h3DVtxZetaPhi", "", kTH3D, {{20, -10, 10}, {16, -0.8, +0.8}, {100, 0., o2::constants::math::TwoPI}});

    // 2D histograms of nSigma
    // before cut
    histos.add("h2DnsigmaPionTpcVsPtBeforeCut", "2D hist of nSigmaTPC vs. pT (pion)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaKaonTpcVsPtBeforeCut", "2D hist of nSigmaTPC vs. pT (kaon)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaProtonTpcVsPtBeforeCut", "2D hist of nSigmaTPC vs. pT (proton)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaPionTofVsPtBeforeCut", "2D hist of nSigmaTOF vs. pT (pion)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaKaonTofVsPtBeforeCut", "2D hist of nSigmaTOF vs. pT (kaon)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaProtonTofVsPtBeforeCut", "2D hist of nSigmaTOF vs. pT (proton)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaPionTpcVsTofBeforeCut", "2D hist of nSigmaTPC vs. nSigmaTOF (pion)", kTH2F, {nSigmaAxis, nSigmaAxis});
    histos.add("h2DnsigmaKaonTpcVsTofBeforeCut", "2D hist of nSigmaTPC vs. nSigmaTOF (kaon)", kTH2F, {nSigmaAxis, nSigmaAxis});
    histos.add("h2DnsigmaProtonTpcVsTofBeforeCut", "2D hist of nSigmaTPC vs. nSigmaTOF (proton)", kTH2F, {nSigmaAxis, nSigmaAxis});
    // after cut
    histos.add("h2DnsigmaPionTpcVsPtAfterCut", "2D hist of nSigmaTPC vs. pT (pion)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaKaonTpcVsPtAfterCut", "2D hist of nSigmaTPC vs. pT (kaon)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaProtonTpcVsPtAfterCut", "2D hist of nSigmaTPC vs. pT (proton)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaPionTofVsPtAfterCut", "2D hist of nSigmaTOF vs. pT (pion)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaKaonTofVsPtAfterCut", "2D hist of nSigmaTOF vs. pT (kaon)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaProtonTofVsPtAfterCut", "2D hist of nSigmaTOF vs. pT (proton)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaPionTpcVsTofAfterCut", "2D hist of nSigmaTPC vs. nSigmaTOF (pion)", kTH2F, {nSigmaAxis, nSigmaAxis});
    histos.add("h2DnsigmaKaonTpcVsTofAfterCut", "2D hist of nSigmaTPC vs. nSigmaTOF (kaon)", kTH2F, {nSigmaAxis, nSigmaAxis});
    histos.add("h2DnsigmaProtonTpcVsTofAfterCut", "2D hist of nSigmaTPC vs. nSigmaTOF (proton)", kTH2F, {nSigmaAxis, nSigmaAxis});

    // Analysis profiles

    histos.add("Prof_A_had", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_C_had", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_D_had", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Bone_had", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Btwo_had", "", {HistType::kTProfile2D, {centAxis, noAxis}});

    histos.add("Prof_A_pi", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_C_pi", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_D_pi", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Bone_pi", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Btwo_pi", "", {HistType::kTProfile2D, {centAxis, noAxis}});

    histos.add("Prof_A_ka", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_C_ka", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_D_ka", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Bone_ka", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Btwo_ka", "", {HistType::kTProfile2D, {centAxis, noAxis}});

    histos.add("Prof_A_prot", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_C_prot", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_D_prot", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Bone_prot", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Btwo_prot", "", {HistType::kTProfile2D, {centAxis, noAxis}});

    // Analysis profile for v02(pT)
    histos.add("Prof_XY", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_XYZ_had", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_Z_had", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_XYZ_pi", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_Z_pi", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_XYZ_ka", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_Z_ka", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_XYZ_prot", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_Z_prot", "", {HistType::kTProfile2D, {centAxis, ptAxis}});

    // check with different normalization for event averaging
    if (cfgV02WeightedFill) {
      histos.add("Prof_XY_weighted", "", {HistType::kTProfile2D, {centAxis, noAxis}});
      histos.add("Prof_XYZ_weighted_had", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
      histos.add("Prof_Z_weighted_had", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    }

    // initial array
    subSample.resize(cfgNSubsample);
    subSampleV02.resize(cfgNSubsample);
    subSampleV02_weighted.resize(cfgNSubsample);
    for (int i = 0; i < cfgNSubsample; i++) {
      subSample[i].resize(20);
      subSampleV02[i].resize(9);
      subSampleV02_weighted[i].resize(9);
    }
    for (int i = 0; i < cfgNSubsample; i++) {
      subSample[i][0] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_A_had", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][1] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_C_had", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][2] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_D_had", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][3] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Bone_had", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][4] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Btwo_had", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));

      subSample[i][5] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_A_pi", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][6] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_C_pi", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][7] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_D_pi", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][8] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Bone_pi", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][9] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Btwo_pi", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));

      subSample[i][10] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_A_ka", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][11] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_C_ka", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][12] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_D_ka", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][13] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Bone_ka", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][14] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Btwo_ka", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));

      subSample[i][15] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_A_prot", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][16] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_C_prot", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][17] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_D_prot", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][18] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Bone_prot", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][19] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Btwo_prot", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));

      subSampleV02[i][0] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_%d/Prof_XY", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSampleV02[i][1] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_%d/Prof_XYZ_had", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSampleV02[i][2] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_%d/Prof_Z_had", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSampleV02[i][3] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_%d/Prof_XYZ_pi", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSampleV02[i][4] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_%d/Prof_Z_pi", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSampleV02[i][5] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_%d/Prof_XYZ_ka", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSampleV02[i][6] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_%d/Prof_Z_ka", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSampleV02[i][7] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_%d/Prof_XYZ_prot", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSampleV02[i][8] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_%d/Prof_Z_prot", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));

      if (cfgV02WeightedFill) {
        subSampleV02_weighted[i][0] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_weighted_%d/Prof_XY_weighted", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
        subSampleV02_weighted[i][1] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_weighted_%d/Prof_XYZ_weighted_had", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
        subSampleV02_weighted[i][2] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_weighted_%d/Prof_Z_weighted_had", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
        subSampleV02_weighted[i][3] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_weighted_%d/Prof_XYZ_weighted_pi", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
        subSampleV02_weighted[i][4] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_weighted_%d/Prof_Z_weighted_pi", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
        subSampleV02_weighted[i][5] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_weighted_%d/Prof_XYZ_weighted_ka", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
        subSampleV02_weighted[i][6] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_weighted_%d/Prof_Z_weighted_ka", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
        subSampleV02_weighted[i][7] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_weighted_%d/Prof_XYZ_weighted_prot", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
        subSampleV02_weighted[i][8] = std::get<std::shared_ptr<TProfile2D>>(histosAnalysis.add(Form("subSampleV02_weighted_%d/Prof_Z_weighted_prot", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      }
    }

    if (cfgEvSelMultCorrelation) {
      cfgFuncParas.multT0CCutPars = cfgFuncParas.cfgMultT0CCutPars;
      cfgFuncParas.multPVT0CCutPars = cfgFuncParas.cfgMultPVT0CCutPars;
      cfgFuncParas.multGlobalPVCutPars = cfgFuncParas.cfgMultGlobalPVCutPars;
      cfgFuncParas.multMultV0ACutPars = cfgFuncParas.cfgMultMultV0ACutPars;
      cfgFuncParas.fMultPVT0CCutLow = new TF1("fMultPVT0CCutLow", cfgFuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultPVT0CCutLow->SetParameters(&(cfgFuncParas.multPVT0CCutPars[0]));
      cfgFuncParas.fMultPVT0CCutHigh = new TF1("fMultPVT0CCutHigh", cfgFuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultPVT0CCutHigh->SetParameters(&(cfgFuncParas.multPVT0CCutPars[0]));
      cfgFuncParas.fMultT0CCutLow = new TF1("fMultT0CCutLow", cfgFuncParas.cfgMultCentLowCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultT0CCutLow->SetParameters(&(cfgFuncParas.multT0CCutPars[0]));
      cfgFuncParas.fMultT0CCutHigh = new TF1("fMultT0CCutHigh", cfgFuncParas.cfgMultCentHighCutFunction->c_str(), 0, 100);
      cfgFuncParas.fMultT0CCutHigh->SetParameters(&(cfgFuncParas.multT0CCutPars[0]));
      cfgFuncParas.fMultGlobalPVCutLow = new TF1("fMultGlobalPVCutLow", cfgFuncParas.cfgMultMultPVLowCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultGlobalPVCutLow->SetParameters(&(cfgFuncParas.multGlobalPVCutPars[0]));
      cfgFuncParas.fMultGlobalPVCutHigh = new TF1("fMultGlobalPVCutHigh", cfgFuncParas.cfgMultMultPVHighCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultGlobalPVCutHigh->SetParameters(&(cfgFuncParas.multGlobalPVCutPars[0]));
      cfgFuncParas.fMultMultV0ACutLow = new TF1("fMultMultV0ACutLow", cfgFuncParas.cfgMultMultV0ALowCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultMultV0ACutLow->SetParameters(&(cfgFuncParas.multMultV0ACutPars[0]));
      cfgFuncParas.fMultMultV0ACutHigh = new TF1("fMultMultV0ACutHigh", cfgFuncParas.cfgMultMultV0AHighCutFunction->c_str(), 0, 4000);
      cfgFuncParas.fMultMultV0ACutHigh->SetParameters(&(cfgFuncParas.multMultV0ACutPars[0]));
      cfgFuncParas.fT0AV0AMean = new TF1("fT0AV0AMean", "[0]+[1]*x", 0, 200000);
      cfgFuncParas.fT0AV0AMean->SetParameters(-1601.0581, 9.417652e-01);
      cfgFuncParas.fT0AV0ASigma = new TF1("fT0AV0ASigma", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 200000);
      cfgFuncParas.fT0AV0ASigma->SetParameters(463.4144, 6.796509e-02, -9.097136e-07, 7.971088e-12, -2.600581e-17);
    }

    if (cfgUsePtDepDCAxy) {
      fPtDepDCAxy = new TF1("ptDepDCAxy", Form("%s", cfgDCAxyFunc->c_str()), 0.001, 1000);
    }
    if (cfgUsePtDepDCAz) {
      fPtDepDCAz = new TF1("ptDepDCAz", Form("%s", cfgDCAzFunc->c_str()), 0.001, 1000);
    }

    if (cfgIsMC) {
      // MC event counts
      histos.add("MCGenerated/hMC", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      histos.add("MCGenerated/hPtEtaPhiCharged", "MC charged particles' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});
      histos.add("MCGenerated/hPtEtaPhiPion", "MC charged pions' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});
      histos.add("MCGenerated/hPtEtaPhiKaon", "MC charged kaons' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});
      histos.add("MCGenerated/hPtEtaPhiProton", "MC charged protons' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});

      histos.add("MCReconstructed/hPtEtaPhiChargedParticle", "MC reconstructed charged particles' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});
      histos.add("MCReconstructed/hPtEtaPhiChargedTrack", "MC reconstructed charged tracks' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});

      histos.add("MCReconstructed/hPtEtaPhiPionParticle", "MC reconstructed pion particles' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});
      histos.add("MCReconstructed/hPtEtaPhiPionTrack", "MC reconstructed pion tracks' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});
      histos.add("MCReconstructed/hPtEtaPhiTruePionTrack", "MC reconstructed pdgcode matched pion tracks' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});

      histos.add("MCReconstructed/hPtEtaPhiKaonParticle", "MC reconstructed kaon particles' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});
      histos.add("MCReconstructed/hPtEtaPhiKaonTrack", "MC reconstructed kaon tracks' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});
      histos.add("MCReconstructed/hPtEtaPhiTrueKaonTrack", "MC reconstructed pdgcode matched kaon tracks' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});

      histos.add("MCReconstructed/hPtEtaPhiProtonParticle", "MC reconstructed proton particles' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});
      histos.add("MCReconstructed/hPtEtaPhiProtonTrack", "MC reconstructed proton tracks' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});
      histos.add("MCReconstructed/hPtEtaPhiTrueProtonTrack", "MC reconstructed pdgcode matched proton tracks' pt, eta, phi", kTH3D, {ptAxis, {100, 0., o2::constants::math::TwoPI}, {100, -2.01, 2.01}});
    }

  } // end init

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  template <typename T>
  bool selectionProton(const T& candidate)
  {
    if (!candidate.hasTPC())
      return false;
    int flag = 0; //! pid check main flag

    if (candidate.pt() > cfgCutPtLower && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC) {
        flag = 1;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < cfgnSigmaCutTOF) {
        flag = 1;
      }
    }
    if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < cfgCutPtUpperPID) {
      const float combNSigmaPr = std::sqrt(std::pow(candidate.tpcNSigmaPr(), 2.0) + std::pow(candidate.tofNSigmaPr(), 2.0));
      const float combNSigmaPi = std::sqrt(std::pow(candidate.tpcNSigmaPi(), 2.0) + std::pow(candidate.tofNSigmaPi(), 2.0));
      const float combNSigmaKa = std::sqrt(std::pow(candidate.tpcNSigmaKa(), 2.0) + std::pow(candidate.tofNSigmaKa(), 2.0));

      int flag2 = 0;
      if (combNSigmaPr < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaPi < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaKa < cfgnSigmaOtherParticles)
        flag2 += 1;

      if (cfgUseNewSeperationPid) {
        if (std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPCHigherPt && std::abs(candidate.tofNSigmaPr()) < cfgnSigmaCutTOFHigherPt) {
          if (!(flag2 > 1) && std::abs(candidate.tpcNSigmaPi()) > cfgnSigmaSeperationCut && std::abs(candidate.tofNSigmaPi()) > cfgnSigmaSeperationCut && std::abs(candidate.tpcNSigmaKa()) > cfgnSigmaSeperationCut && std::abs(candidate.tofNSigmaKa()) > cfgnSigmaSeperationCut)
            flag = 1;
        }
      } else {
        if (!(flag2 > 1) && !(combNSigmaPr > combNSigmaPi) && !(combNSigmaPr > combNSigmaKa)) {
          if (combNSigmaPr < cfgnSigmaCutCombTPCTOF) {
            flag = 1;
          }
        }
      }
    }
    if (flag == 1)
      return true;
    else
      return false;
  }

  template <typename T>
  bool selectionPion(const T& candidate)
  {
    if (!candidate.hasTPC())
      return false;
    int flag = 0; //! pid check main flag

    if (candidate.pt() > cfgCutPtLower && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPC) {
        flag = 1;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaPi()) < cfgnSigmaCutTOF) {
        flag = 1;
      }
    }
    if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < cfgCutPtUpperPID) {
      const float combNSigmaPr = std::sqrt(std::pow(candidate.tpcNSigmaPr(), 2.0) + std::pow(candidate.tofNSigmaPr(), 2.0));
      const float combNSigmaPi = std::sqrt(std::pow(candidate.tpcNSigmaPi(), 2.0) + std::pow(candidate.tofNSigmaPi(), 2.0));
      const float combNSigmaKa = std::sqrt(std::pow(candidate.tpcNSigmaKa(), 2.0) + std::pow(candidate.tofNSigmaKa(), 2.0));

      int flag2 = 0;
      if (combNSigmaPr < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaPi < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaKa < cfgnSigmaOtherParticles)
        flag2 += 1;

      if (cfgUseNewSeperationPid) {
        if (std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPCHigherPt && std::abs(candidate.tofNSigmaPi()) < cfgnSigmaCutTOFHigherPt) {
          if (!(flag2 > 1) && std::abs(candidate.tpcNSigmaKa()) > cfgnSigmaSeperationCut && std::abs(candidate.tofNSigmaKa()) > cfgnSigmaSeperationCut && std::abs(candidate.tpcNSigmaPr()) > cfgnSigmaSeperationCut && std::abs(candidate.tofNSigmaPr()) > cfgnSigmaSeperationCut)
            flag = 1;
        }
      } else {
        if (!(flag2 > 1) && !(combNSigmaPi > combNSigmaPr) && !(combNSigmaPi > combNSigmaKa)) {
          if (combNSigmaPi < cfgnSigmaCutCombTPCTOF) {
            flag = 1;
          }
        }
      }
    }
    if (flag == 1)
      return true;
    else
      return false;
  }

  template <typename T>
  bool selectionKaon(const T& candidate)
  {
    if (!candidate.hasTPC())
      return false;
    int flag = 0; //! pid check main flag

    if (candidate.pt() > cfgCutPtLower && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPC) {
        flag = 1;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaKa()) < cfgnSigmaCutTOF) {
        flag = 1;
      }
    }
    if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < cfgCutPtUpperPID) {
      const float combNSigmaPr = std::sqrt(std::pow(candidate.tpcNSigmaPr(), 2.0) + std::pow(candidate.tofNSigmaPr(), 2.0));
      const float combNSigmaPi = std::sqrt(std::pow(candidate.tpcNSigmaPi(), 2.0) + std::pow(candidate.tofNSigmaPi(), 2.0));
      const float combNSigmaKa = std::sqrt(std::pow(candidate.tpcNSigmaKa(), 2.0) + std::pow(candidate.tofNSigmaKa(), 2.0));

      int flag2 = 0;
      if (combNSigmaPr < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaPi < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaKa < cfgnSigmaOtherParticles)
        flag2 += 1;

      if (cfgUseNewSeperationPid) {
        if (std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPCHigherPt && std::abs(candidate.tofNSigmaKa()) < cfgnSigmaCutTOFHigherPt) {
          if (!(flag2 > 1) && std::abs(candidate.tpcNSigmaPi()) > cfgnSigmaSeperationCut && std::abs(candidate.tofNSigmaPi()) > cfgnSigmaSeperationCut && std::abs(candidate.tpcNSigmaPr()) > cfgnSigmaSeperationCut && std::abs(candidate.tofNSigmaPr()) > cfgnSigmaSeperationCut)
            flag = 1;
        }
      } else {
        if (!(flag2 > 1) && !(combNSigmaKa > combNSigmaPi) && !(combNSigmaKa > combNSigmaPr)) {
          if (combNSigmaKa < cfgnSigmaCutCombTPCTOF) {
            flag = 1;
          }
        }
      }
    }
    if (flag == 1)
      return true;
    else
      return false;
  }

  template <typename TTrack>
  int getNsigmaPID(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {itsResponse.nSigmaITS<o2::track::PID::Pion>(track), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = 0; // 0 = not identified, 1 = pion, 2 = kaon, 3 = proton

    std::array<float, 3> nSigmaToUse = cfgUseItsPID ? nSigmaITS : nSigmaTPC;             // Choose which nSigma to use: TPC or ITS
    std::array<float, 6> detectorNsigmaCut = cfgUseItsPID ? itsNsigmaCut : tpcNsigmaCut; // Choose which nSigma to use: TPC or ITS

    bool isPion, isKaon, isProton;
    bool isDetectedPion = nSigmaToUse[PIONS] < detectorNsigmaCut[kPionUpCut] && nSigmaToUse[PIONS] > detectorNsigmaCut[kPionLowCut];
    bool isDetectedKaon = nSigmaToUse[KAONS] < detectorNsigmaCut[kKaonUpCut] && nSigmaToUse[KAONS] > detectorNsigmaCut[kKaonLowCut];
    bool isDetectedProton = nSigmaToUse[PROTONS] < detectorNsigmaCut[kProtonUpCut] && nSigmaToUse[PROTONS] > detectorNsigmaCut[kProtonLowCut];

    bool isTofPion = nSigmaTOF[PIONS] < tofNsigmaCut[kPionUpCut] && nSigmaTOF[PIONS] > tofNsigmaCut[kPionLowCut];
    bool isTofKaon = nSigmaTOF[KAONS] < tofNsigmaCut[kKaonUpCut] && nSigmaTOF[KAONS] > tofNsigmaCut[kKaonLowCut];
    bool isTofProton = nSigmaTOF[PROTONS] < tofNsigmaCut[kProtonUpCut] && nSigmaTOF[PROTONS] > tofNsigmaCut[kProtonLowCut];

    if (track.pt() > cfgPtCutTOF && !track.hasTOF()) {
      return 0;
    } else if (track.pt() > cfgPtCutTOF && track.hasTOF()) {
      isPion = isTofPion && isDetectedPion;
      isKaon = isTofKaon && isDetectedKaon;
      isProton = isTofProton && isDetectedProton;
    } else {
      isPion = isDetectedPion;
      isKaon = isDetectedKaon;
      isProton = isDetectedProton;
    }

    if ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton)) {
      return 0; // more than one particle satisfy the criteria
    }

    if (isPion) {
      pid = PIONS + 1;
    } else if (isKaon) {
      pid = KAONS + 1;
    } else if (isProton) {
      pid = PROTONS + 1;
    } else {
      return 0; // no particle satisfies the criteria
    }

    return pid; // 0 = not identified, 1 = pion, 2 = kaon, 3 = proton
  }

  // additional multiplicity correlation based event selection cuts
  template <typename TCollision>
  bool eventSelectedSmallion(TCollision collision, const int multTrk, const float centrality)
  {
    auto multNTracksPV = collision.multNTracksPV();

    if (cfgEvSelMultCorrelation) {
      if (cfgFuncParas.cfgMultPVT0CCutEnabled) {
        if (multNTracksPV < cfgFuncParas.fMultPVT0CCutLow->Eval(centrality))
          return 0;
        if (multNTracksPV > cfgFuncParas.fMultPVT0CCutHigh->Eval(centrality))
          return 0;
      }

      if (cfgFuncParas.cfgMultT0CCutEnabled) {
        if (multTrk < cfgFuncParas.fMultT0CCutLow->Eval(centrality))
          return 0;
        if (multTrk > cfgFuncParas.fMultT0CCutHigh->Eval(centrality))
          return 0;
      }

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
    }

    float sigma = 5.0;
    if (cfgEvSelV0AT0ACut && (std::fabs(collision.multFV0A() - cfgFuncParas.fT0AV0AMean->Eval(collision.multFT0A())) > sigma * cfgFuncParas.fT0AV0ASigma->Eval(collision.multFT0A())))
      return 0;

    return 1;
  }

  template <typename TCollision>
  bool eventSelectionDefaultCuts(TCollision coll)
  {
    histos.fill(HIST("hEventStatData"), 0.5);
    if (!coll.sel8()) {
      return 0;
    }

    histos.fill(HIST("hEventStatData"), 1.5);
    if (cfgUseGoodITSLayerAllCut && !(coll.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))) {
      return 0;
    }

    histos.fill(HIST("hEventStatData"), 2.5);
    if (cfgEvSelkNoSameBunchPileup && !(coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
      return 0;
    }

    histos.fill(HIST("hEventStatData"), 3.5);
    if (cfgEvSelkNoITSROFrameBorder && !(coll.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
      return 0;
    }

    histos.fill(HIST("hEventStatData"), 4.5);
    if (cfgEvSelkNoTimeFrameBorder && !(coll.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))) {
      return 0;
    }

    histos.fill(HIST("hEventStatData"), 5.5);
    if (cfgEvSelUseGoodZvtxFT0vsPV && !(coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return 0;
    }

    histos.fill(HIST("hEventStatData"), 6.5);
    // events with selection bits based on occupancy time pattern
    if (cfgEvSelUseOcuppancyTimeCut && !(coll.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
      return 0;
    }

    histos.fill(HIST("hEventStatData"), 7.5);
    int occupancy = coll.trackOccupancyInTimeRange();
    if (cfgEvSelSetOcuppancyRange && (occupancy < cfgMinOccupancy || occupancy > cfgMaxOccupancy)) {
      return 0;
    }

    histos.fill(HIST("hEventStatData"), 8.5);
    return 1;
  }

  template <typename C, typename T>
  void fillMultCorrPlotsBeforeSel(C const& coll, T const& inputTracks)
  {
    histos.fill(HIST("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_PVTracks_beforeSel"), coll.multNTracksPV(), inputTracks.size());
    histos.fill(HIST("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_centFT0C_beforeSel"), coll.centFT0C(), inputTracks.size());
    histos.fill(HIST("MultCorrelationPlots/BeforeSelection/His2D_PVTracks_centFT0C_beforeSel"), coll.centFT0C(), coll.multNTracksPV());
    histos.fill(HIST("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_V0ATracks_beforeSel"), coll.multFV0A(), inputTracks.size());
    histos.fill(HIST("MultCorrelationPlots/BeforeSelection/His2D_globalTracks_T0ATracks_beforeSel"), coll.multFT0A(), inputTracks.size());
    histos.fill(HIST("MultCorrelationPlots/BeforeSelection/His2D_V0ATracks_T0CTracks_beforeSel"), coll.multFT0C(), coll.multFV0A());
  }

  template <typename C, typename T>
  void fillMultCorrPlotsAfterSel(C const& coll, T const& inputTracks)
  {
    histos.fill(HIST("MultCorrelationPlots/AfterSelection/His2D_globalTracks_PVTracks_afterSel"), coll.multNTracksPV(), inputTracks.size());
    histos.fill(HIST("MultCorrelationPlots/AfterSelection/His2D_globalTracks_centFT0C_afterSel"), coll.centFT0C(), inputTracks.size());
    histos.fill(HIST("MultCorrelationPlots/AfterSelection/His2D_PVTracks_centFT0C_afterSel"), coll.centFT0C(), coll.multNTracksPV());
    histos.fill(HIST("MultCorrelationPlots/AfterSelection/His2D_globalTracks_V0ATracks_afterSel"), coll.multFV0A(), inputTracks.size());
    histos.fill(HIST("MultCorrelationPlots/AfterSelection/His2D_globalTracks_T0ATracks_afterSel"), coll.multFT0A(), inputTracks.size());
    histos.fill(HIST("MultCorrelationPlots/AfterSelection/His2D_V0ATracks_T0CTracks_afterSel"), coll.multFT0C(), coll.multFV0A());
  }

  template <typename T>
  float getPhiWeight(const T& candidate, float vtxz)
  {
    if (!cfgLoadPhiWeights || !hWeightPhiFunctionVzEtaPhi) {
      return 1.0;
    }
    int bin = hWeightPhiFunctionVzEtaPhi->FindBin(vtxz, candidate.eta(), candidate.phi());
    float weight = hWeightPhiFunctionVzEtaPhi->GetBinContent(bin);
    if (!std::isfinite(weight) || weight <= 0) {
      return 1.0;
    }
    return weight;
  }

  template <typename T>
  float getEffAllCharged(const T& candidate)
  {
    if (!cfgLoadPtEffWeights || !hEffAllCharged) {
      return 1.0;
    }
    int bin = hEffAllCharged->FindBin(candidate.pt());
    float eff = hEffAllCharged->GetBinContent(bin);
    float ptweight = 1.0 / eff;
    if (!std::isfinite(ptweight) || ptweight <= 0) {
      return 1.0;
    }
    return ptweight;
  }

  template <typename T>
  float getEffPion(const T& candidate)
  {
    if (!cfgLoadPtEffWeights || !hEffPion) {
      return 1.0;
    }
    int bin = hEffPion->FindBin(candidate.pt());
    float eff = hEffPion->GetBinContent(bin);
    float ptweight = 1.0 / eff;
    if (!std::isfinite(ptweight) || ptweight <= 0) {
      return 1.0;
    }
    return ptweight;
  }

  template <typename T>
  float getEffKaon(const T& candidate)
  {
    if (!cfgLoadPtEffWeights || !hEffKaon) {
      return 1.0;
    }
    int bin = hEffKaon->FindBin(candidate.pt());
    float eff = hEffKaon->GetBinContent(bin);
    float ptweight = 1.0 / eff;
    if (!std::isfinite(ptweight) || ptweight <= 0) {
      return 1.0;
    }
    return ptweight;
  }

  template <typename T>
  float getEffProton(const T& candidate)
  {
    if (!cfgLoadPtEffWeights || !hEffProton) {
      return 1.0;
    }
    int bin = hEffProton->FindBin(candidate.pt());
    float eff = hEffProton->GetBinContent(bin);
    float ptweight = 1.0 / eff;
    if (!std::isfinite(ptweight) || ptweight <= 0) {
      return 1.0;
    }
    return ptweight;
  }

  // process MC recosnstructed data
  void processMCRec(MyMCRecCollisions::iterator const& collision, MyMCTracks const& tracks, aod::McCollisions const&, aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("MCGenerated/hMC"), 5.5);

    if (!collision.has_mcCollision()) {
      return;
    }
    histos.fill(HIST("MCGenerated/hMC"), 6.5);

    if (!eventSelectionDefaultCuts(collision)) {
      return;
    }
    histos.fill(HIST("MCGenerated/hMC"), 7.5);

    fillMultCorrPlotsBeforeSel(collision, tracks);

    const auto centralityFT0C = collision.centFT0C();
    if (cfgUseSmallIonAdditionalEventCut && !eventSelectedSmallion(collision, tracks.size(), centralityFT0C))
      return;

    if (cfgUseSmallIonAdditionalEventCut) {
      fillMultCorrPlotsAfterSel(collision, tracks);
    }

    histos.fill(HIST("MCGenerated/hMC"), 8.5);

    // Centrality
    double cent = 0.0;
    if (cfgCentralityChoice == kFT0C)
      cent = collision.centFT0C();
    else if (cfgCentralityChoice == kFT0A)
      cent = collision.centFT0A();
    else if (cfgCentralityChoice == kFT0M)
      cent = collision.centFT0M();
    else if (cfgCentralityChoice == kFV0A)
      cent = collision.centFV0A();

    histos.fill(HIST("hZvtx_after_sel"), collision.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), collision.multNTracksPV(), tracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), tracks.size(), centralityFT0C);

    // Calculating generated level pt distributions
    auto mcColl = collision.mcCollision();
    // Slice particles belonging only to this MC collision
    auto particlesThisEvent = mcParticles.sliceBy(perMcCollision, mcColl.globalIndex());

    for (const auto& mcParticle : particlesThisEvent) {
      if (!mcParticle.has_mcCollision())
        continue;

      // charged check
      auto pdgEntry = TDatabasePDG::Instance()->GetParticle(mcParticle.pdgCode());
      if (!pdgEntry)
        continue;
      if (pdgEntry->Charge() == 0)
        continue;

      if (mcParticle.isPhysicalPrimary()) {
        if ((mcParticle.pt() > cfgCutPtLower) && (mcParticle.pt() < cfgCutPtUpper) && (std::abs(mcParticle.eta()) < cfgCutEta)) {
          histos.fill(HIST("MCGenerated/hPtEtaPhiCharged"), mcParticle.pt(), mcParticle.eta(), mcParticle.phi());

          auto pdgcode = std::abs(mcParticle.pdgCode());

          if (pdgcode == PDG_t::kPiPlus)
            histos.fill(HIST("MCGenerated/hPtEtaPhiPion"), mcParticle.pt(), mcParticle.eta(), mcParticle.phi());

          if (pdgcode == PDG_t::kKPlus)
            histos.fill(HIST("MCGenerated/hPtEtaPhiKaon"), mcParticle.pt(), mcParticle.eta(), mcParticle.phi());

          if (pdgcode == PDG_t::kProton)
            histos.fill(HIST("MCGenerated/hPtEtaPhiProton"), mcParticle.pt(), mcParticle.eta(), mcParticle.phi());
        }
      }
    } //! end mc particle loop

    for (const auto& track : tracks) { // Loop over tracks

      if (!track.has_collision()) {
        continue;
      }
      if (!track.has_mcParticle()) { //! check if track has corresponding MC particle
        continue;
      }

      auto particle = track.mcParticle();
      if (!particle.has_mcCollision()) {
        continue;
      }

      if (!track.isPVContributor()) {
        continue;
      }

      if (!(track.itsNCls() > cfgITScluster) || !(track.tpcNClsFound() >= cfgTPCcluster) || !(track.tpcNClsCrossedRows() >= cfgTPCnCrossedRows)) {
        continue;
      }

      if (cfgUsePtDepDCAxy && !(std::abs(track.dcaXY()) < fPtDepDCAxy->Eval(track.pt()))) {
        continue;
      }
      if (cfgUsePtDepDCAz && !(std::abs(track.dcaZ()) < fPtDepDCAz->Eval(track.pt()))) {
        continue;
      }

      if (track.sign() == 0)
        continue;

      if (particle.isPhysicalPrimary()) {
        if ((particle.pt() > cfgCutPtLower) && (particle.pt() < cfgCutPtUpper) && (std::abs(particle.eta()) < cfgCutEta)) {

          histos.fill(HIST("MCReconstructed/hPtEtaPhiChargedParticle"), particle.pt(), particle.eta(), particle.phi());
          histos.fill(HIST("MCReconstructed/hPtEtaPhiChargedTrack"), track.pt(), track.eta(), track.phi());

          // PID QAs before selection
          double nSigmaTpcPi = track.tpcNSigmaPi();
          double nSigmaTpcKa = track.tpcNSigmaKa();
          double nSigmaTpcProt = track.tpcNSigmaPr();
          double nSigmaTofPi = track.tofNSigmaPi();
          double nSigmaTofKa = track.tofNSigmaKa();
          double nSigmaTofProt = track.tofNSigmaPr();
          histos.fill(HIST("h2DnsigmaPionTpcVsPtBeforeCut"), track.pt(), nSigmaTpcPi);
          histos.fill(HIST("h2DnsigmaKaonTpcVsPtBeforeCut"), track.pt(), nSigmaTpcKa);
          histos.fill(HIST("h2DnsigmaProtonTpcVsPtBeforeCut"), track.pt(), nSigmaTpcProt);
          histos.fill(HIST("h2DnsigmaPionTofVsPtBeforeCut"), track.pt(), nSigmaTofPi);
          histos.fill(HIST("h2DnsigmaKaonTofVsPtBeforeCut"), track.pt(), nSigmaTofKa);
          histos.fill(HIST("h2DnsigmaProtonTofVsPtBeforeCut"), track.pt(), nSigmaTofProt);
          histos.fill(HIST("h2DnsigmaPionTpcVsTofBeforeCut"), nSigmaTpcPi, nSigmaTofPi);
          histos.fill(HIST("h2DnsigmaKaonTpcVsTofBeforeCut"), nSigmaTpcKa, nSigmaTofKa);
          histos.fill(HIST("h2DnsigmaProtonTpcVsTofBeforeCut"), nSigmaTpcProt, nSigmaTofProt);

          // identified particles selection
          bool isPion = false;
          bool isKaon = false;
          bool isProton = false;

          if (cfgUseRun3V2PID) {
            int pidVal = getNsigmaPID(track);
            if (pidVal == PIONS + 1)
              isPion = true;
            if (pidVal == KAONS + 1)
              isKaon = true;
            if (pidVal == PROTONS + 1)
              isProton = true;
          } else {
            isPion = selectionPion(track);
            isKaon = selectionKaon(track);
            isProton = selectionProton(track);
          }

          // PID QAs after selection
          if (isPion) {
            histos.fill(HIST("h2DnsigmaPionTpcVsPtAfterCut"), track.pt(), nSigmaTpcPi);
            histos.fill(HIST("h2DnsigmaPionTofVsPtAfterCut"), track.pt(), nSigmaTofPi);
            histos.fill(HIST("h2DnsigmaPionTpcVsTofAfterCut"), nSigmaTpcPi, nSigmaTofPi);
          }
          if (isKaon) {
            histos.fill(HIST("h2DnsigmaKaonTpcVsPtAfterCut"), track.pt(), nSigmaTpcKa);
            histos.fill(HIST("h2DnsigmaKaonTofVsPtAfterCut"), track.pt(), nSigmaTofKa);
            histos.fill(HIST("h2DnsigmaKaonTpcVsTofAfterCut"), nSigmaTpcKa, nSigmaTofKa);
          }
          if (isProton) {
            histos.fill(HIST("h2DnsigmaProtonTpcVsPtAfterCut"), track.pt(), nSigmaTpcProt);
            histos.fill(HIST("h2DnsigmaProtonTofVsPtAfterCut"), track.pt(), nSigmaTofProt);
            histos.fill(HIST("h2DnsigmaProtonTpcVsTofAfterCut"), nSigmaTpcProt, nSigmaTofProt);
          }

          auto pdgcodeRec = std::abs(particle.pdgCode());

          if (isPion) {
            histos.fill(HIST("MCReconstructed/hPtEtaPhiPionParticle"), particle.pt(), particle.eta(), particle.phi());
            histos.fill(HIST("MCReconstructed/hPtEtaPhiPionTrack"), track.pt(), track.eta(), track.phi());
            if (pdgcodeRec == PDG_t::kPiPlus)
              histos.fill(HIST("MCReconstructed/hPtEtaPhiTruePionTrack"), track.pt(), track.eta(), track.phi());
          }

          if (isKaon) {
            histos.fill(HIST("MCReconstructed/hPtEtaPhiKaonParticle"), particle.pt(), particle.eta(), particle.phi());
            histos.fill(HIST("MCReconstructed/hPtEtaPhiKaonTrack"), track.pt(), track.eta(), track.phi());
            if (pdgcodeRec == PDG_t::kKPlus)
              histos.fill(HIST("MCReconstructed/hPtEtaPhiTrueKaonTrack"), track.pt(), track.eta(), track.phi());
          }

          if (isProton) {
            histos.fill(HIST("MCReconstructed/hPtEtaPhiProtonParticle"), particle.pt(), particle.eta(), particle.phi());
            histos.fill(HIST("MCReconstructed/hPtEtaPhiProtonTrack"), track.pt(), track.eta(), track.phi());
            if (pdgcodeRec == PDG_t::kProton)
              histos.fill(HIST("MCReconstructed/hPtEtaPhiTrueProtonTrack"), track.pt(), track.eta(), track.phi());
          }
        }
      }
    } // end track loop
  }
  PROCESS_SWITCH(V0ptHadPiKaProt, processMCRec, "Process Monte-carlo data", false);

  // process Data
  void processData(AodCollisions::iterator const& coll, aod::BCsWithTimestamps const&, AodTracks const& inputTracks)
  {
    if (!eventSelectionDefaultCuts(coll)) {
      return;
    }

    fillMultCorrPlotsBeforeSel(coll, inputTracks);

    const auto centralityFT0C = coll.centFT0C();
    if (cfgUseSmallIonAdditionalEventCut && !eventSelectedSmallion(coll, inputTracks.size(), centralityFT0C))
      return;

    if (cfgUseSmallIonAdditionalEventCut) {
      fillMultCorrPlotsAfterSel(coll, inputTracks);
    }

    // Centrality
    double cent = 0.0;
    if (cfgCentralityChoice == kFT0C)
      cent = coll.centFT0C();
    else if (cfgCentralityChoice == kFT0A)
      cent = coll.centFT0A();
    else if (cfgCentralityChoice == kFT0M)
      cent = coll.centFT0M();
    else if (cfgCentralityChoice == kFV0A)
      cent = coll.centFV0A();

    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), inputTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), inputTracks.size(), cent);

    // Analysis variables for v0(pT)
    int nbinsHad = 20;
    int nbinsPid = 18;
    double binsarray[21] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
    TH1D* fPtProfileHad = new TH1D("fPtProfileHad", "fPtProfileHad", 20, binsarray);
    TH1D* fPtProfilePi = new TH1D("fPtProfilePi", "fPtProfilePi", 20, binsarray);
    TH1D* fPtProfileKa = new TH1D("fPtProfileKa", "fPtProfileKa", 20, binsarray);
    TH1D* fPtProfileProt = new TH1D("fPtProfileProt", "fPtProfileProt", 20, binsarray);
    double pTsumEtaLeftHad = 0.0;
    double nSumEtaLeftHad = 0.0;
    double pTsumEtaRightHad = 0.0;
    double nSumEtaRightHad = 0.0;
    double nSumEtaLeftPi = 0.0;
    double nSumEtaLeftKa = 0.0;
    double nSumEtaLeftProt = 0.0;

    // Analysis variables for v02(pT)
    TH1D* fPtProfileHadInWinB = new TH1D("fPtProfileHadInWinB", "fPtProfileHadInWinB", 20, binsarray);
    TH1D* fPtProfilePiInWinB = new TH1D("fPtProfilePiInWinB", "fPtProfilePiInWinB", 20, binsarray);
    TH1D* fPtProfileKaInWinB = new TH1D("fPtProfileKaInWinB", "fPtProfileKaInWinB", 20, binsarray);
    TH1D* fPtProfileProtInWinB = new TH1D("fPtProfileProtInWinB", "fPtProfileProtInWinB", 20, binsarray);
    double nSumInWinB = 0.0; // for Z = f(pT) = n(pT)/N_B in window B

    double nSumInWinA = 0.0; // for X (in window A) to calculate v2^2
    double nSumInWinC = 0.0; // for Y (in window C) to calculate v2^2
    TComplex vecQInWinA = TComplex(0., 0.);
    TComplex vecQInWinC = TComplex(0., 0.);

    for (const auto& track : inputTracks) { // Loop over tracks

      if (!track.has_collision()) {
        continue;
      }

      if (!track.isPVContributor()) {
        continue;
      }

      if (!(track.itsNCls() > cfgITScluster) || !(track.tpcNClsFound() >= cfgTPCcluster) || !(track.tpcNClsCrossedRows() >= cfgTPCnCrossedRows)) {
        continue;
      }

      if (cfgUsePtDepDCAxy && !(std::abs(track.dcaXY()) < fPtDepDCAxy->Eval(track.pt()))) {
        continue;
      }
      if (cfgUsePtDepDCAz && !(std::abs(track.dcaZ()) < fPtDepDCAz->Eval(track.pt()))) {
        continue;
      }

      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
      histos.fill(HIST("hEta"), track.eta());
      histos.fill(HIST("hPhi"), track.phi());
      histos.fill(HIST("hDcaXY"), track.dcaXY());
      histos.fill(HIST("hDcaZ"), track.dcaZ());

      double trkPt = track.pt();
      double trkEta = track.eta();
      double trkPhi = track.phi();

      // inclusive charged particles
      if (track.sign() != 0) {
        if (trkEta < cfgCutEtaLeft) {
          fPtProfileHad->Fill(trkPt);
          pTsumEtaLeftHad += trkPt;
          nSumEtaLeftHad += 1.0;
        }
        if (trkEta > cfgCutEtaRight) {
          pTsumEtaRightHad += trkPt;
          nSumEtaRightHad += 1.0;
        }
      }

      double phiweight = 1.0;
      if (cfgLoadPhiWeights) {
        phiweight = getPhiWeight(track, coll.posZ()); // NUA weight
      }
      double effweight = 1.0;
      if (cfgLoadPtEffWeights) {
        effweight = 1.0 / getEffAllCharged(track); // NUE weight
      }
      double weight = phiweight * effweight;

      if (track.sign() != 0 && trkPt < cfgCutPtMaxForV02) {
        histos.fill(HIST("h3DVtxZetaPhi"), coll.posZ(), trkEta, trkPhi);

        // fill subevent C for v2^2 in v02(pT)
        if (cfgCutEtaWindowB < trkEta && trkEta < cfgCutEta) {
          vecQInWinC += weight * TComplex(TMath::Cos(2. * trkPhi), TMath::Sin(2. * trkPhi));
          nSumInWinC += weight;
        }

        // fill subevent A for v2^2 in v02(pT)
        if (-1.0 * cfgCutEta < trkEta && trkEta < -1.0 * cfgCutEtaWindowB) {
          vecQInWinA += weight * TComplex(TMath::Cos(2. * trkPhi), TMath::Sin(2. * trkPhi));
          nSumInWinA += weight;
        }

        // fill subevent B for f(pT) in v02(pT)
        if (std::abs(trkEta) < cfgCutEtaWindowB) {
          fPtProfileHadInWinB->Fill(trkPt, effweight);
          nSumInWinB += effweight;
        }
      }

      // PID QAs before selection
      double nSigmaTpcPi = track.tpcNSigmaPi();
      double nSigmaTpcKa = track.tpcNSigmaKa();
      double nSigmaTpcProt = track.tpcNSigmaPr();
      double nSigmaTofPi = track.tofNSigmaPi();
      double nSigmaTofKa = track.tofNSigmaKa();
      double nSigmaTofProt = track.tofNSigmaPr();
      histos.fill(HIST("h2DnsigmaPionTpcVsPtBeforeCut"), trkPt, nSigmaTpcPi);
      histos.fill(HIST("h2DnsigmaKaonTpcVsPtBeforeCut"), trkPt, nSigmaTpcKa);
      histos.fill(HIST("h2DnsigmaProtonTpcVsPtBeforeCut"), trkPt, nSigmaTpcProt);
      histos.fill(HIST("h2DnsigmaPionTofVsPtBeforeCut"), trkPt, nSigmaTofPi);
      histos.fill(HIST("h2DnsigmaKaonTofVsPtBeforeCut"), trkPt, nSigmaTofKa);
      histos.fill(HIST("h2DnsigmaProtonTofVsPtBeforeCut"), trkPt, nSigmaTofProt);
      histos.fill(HIST("h2DnsigmaPionTpcVsTofBeforeCut"), nSigmaTpcPi, nSigmaTofPi);
      histos.fill(HIST("h2DnsigmaKaonTpcVsTofBeforeCut"), nSigmaTpcKa, nSigmaTofKa);
      histos.fill(HIST("h2DnsigmaProtonTpcVsTofBeforeCut"), nSigmaTpcProt, nSigmaTofProt);

      // identified particles selection
      bool isPion = false;
      bool isKaon = false;
      bool isProton = false;

      if (cfgUseRun3V2PID) {
        int pidVal = getNsigmaPID(track);
        if (pidVal == PIONS + 1)
          isPion = true;
        if (pidVal == KAONS + 1)
          isKaon = true;
        if (pidVal == PROTONS + 1)
          isProton = true;
      } else {
        isPion = selectionPion(track);
        isKaon = selectionKaon(track);
        isProton = selectionProton(track);
      }

      // PID QAs after selection
      if (isPion) {
        histos.fill(HIST("h2DnsigmaPionTpcVsPtAfterCut"), trkPt, nSigmaTpcPi);
        histos.fill(HIST("h2DnsigmaPionTofVsPtAfterCut"), trkPt, nSigmaTofPi);
        histos.fill(HIST("h2DnsigmaPionTpcVsTofAfterCut"), nSigmaTpcPi, nSigmaTofPi);
      }
      if (isKaon) {
        histos.fill(HIST("h2DnsigmaKaonTpcVsPtAfterCut"), trkPt, nSigmaTpcKa);
        histos.fill(HIST("h2DnsigmaKaonTofVsPtAfterCut"), trkPt, nSigmaTofKa);
        histos.fill(HIST("h2DnsigmaKaonTpcVsTofAfterCut"), nSigmaTpcKa, nSigmaTofKa);
      }
      if (isProton) {
        histos.fill(HIST("h2DnsigmaProtonTpcVsPtAfterCut"), trkPt, nSigmaTpcProt);
        histos.fill(HIST("h2DnsigmaProtonTofVsPtAfterCut"), trkPt, nSigmaTofProt);
        histos.fill(HIST("h2DnsigmaProtonTpcVsTofAfterCut"), nSigmaTpcProt, nSigmaTofProt);
      }

      if (track.sign() != 0) {
        if (trkPt < cfgCutPtUpperPID) {
          if (trkEta < cfgCutEtaLeft) {
            if (isPion) {
              fPtProfilePi->Fill(trkPt);
              nSumEtaLeftPi += 1.0;
            }
            if (isKaon) {
              fPtProfileKa->Fill(trkPt);
              nSumEtaLeftKa += 1.0;
            }
            if (isProton && trkPt > cfgCutPtLowerProt) {
              fPtProfileProt->Fill(trkPt);
              nSumEtaLeftProt += 1.0;
            }
          }
        }
      }

      double effweightPion = 1.0;
      double effweightKaon = 1.0;
      double effweightProton = 1.0;
      if (cfgLoadPtEffWeights) {
        effweightPion = 1.0 / getEffPion(track);     // NUE weight for pion
        effweightKaon = 1.0 / getEffKaon(track);     // NUE weight for kaon
        effweightProton = 1.0 / getEffProton(track); // NUE weight for proton
      }

      // fill subevent B for ***identified particles'*** f(pT) in v02(pT)
      if (track.sign() != 0 && trkPt < cfgCutPtMaxForV02) {
        if (std::abs(trkEta) < cfgCutEtaWindowB) {
          if (isPion) {
            fPtProfilePiInWinB->Fill(trkPt, effweightPion);
          }
          if (isKaon) {
            fPtProfileKaInWinB->Fill(trkPt, effweightKaon);
          }
          if (isProton && trkPt > cfgCutPtLowerProt) {
            fPtProfileProtInWinB->Fill(trkPt, effweightProton);
          }
        }
      }

    } // End track loop

    // selecting subsample and filling profiles
    float lRandom = funRndm->Rndm();
    int sampleIndex = static_cast<int>(cfgNSubsample * lRandom);

    if (nSumEtaRightHad > 0 && nSumEtaLeftHad > 0) {
      for (int i = 0; i < nbinsHad; i++) {
        histos.get<TProfile2D>(HIST("Prof_A_had"))->Fill(cent, fPtProfileHad->GetBinCenter(i + 1), (fPtProfileHad->GetBinContent(i + 1) / nSumEtaLeftHad));
        histos.get<TProfile2D>(HIST("Prof_C_had"))->Fill(cent, fPtProfileHad->GetBinCenter(i + 1), ((fPtProfileHad->GetBinContent(i + 1) / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));
        histos.get<TProfile2D>(HIST("Prof_Bone_had"))->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        histos.get<TProfile2D>(HIST("Prof_Btwo_had"))->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
        histos.get<TProfile2D>(HIST("Prof_D_had"))->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));

        subSample[sampleIndex][0]->Fill(cent, fPtProfileHad->GetBinCenter(i + 1), (fPtProfileHad->GetBinContent(i + 1) / nSumEtaLeftHad));
        subSample[sampleIndex][1]->Fill(cent, fPtProfileHad->GetBinCenter(i + 1), ((fPtProfileHad->GetBinContent(i + 1) / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][2]->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][3]->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        subSample[sampleIndex][4]->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
      }
    }

    if (nSumEtaRightHad > 0 && nSumEtaLeftHad > 0 && nSumEtaLeftPi > 0) {
      for (int i = 0; i < nbinsPid; i++) {
        histos.get<TProfile2D>(HIST("Prof_A_pi"))->Fill(cent, fPtProfilePi->GetBinCenter(i + 1), (fPtProfilePi->GetBinContent(i + 1) / nSumEtaLeftPi));
        histos.get<TProfile2D>(HIST("Prof_C_pi"))->Fill(cent, fPtProfilePi->GetBinCenter(i + 1), ((fPtProfilePi->GetBinContent(i + 1) / nSumEtaLeftPi) * (pTsumEtaRightHad / nSumEtaRightHad)));
        histos.get<TProfile2D>(HIST("Prof_Bone_pi"))->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        histos.get<TProfile2D>(HIST("Prof_Btwo_pi"))->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
        histos.get<TProfile2D>(HIST("Prof_D_pi"))->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));

        subSample[sampleIndex][5]->Fill(cent, fPtProfilePi->GetBinCenter(i + 1), (fPtProfilePi->GetBinContent(i + 1) / nSumEtaLeftPi));
        subSample[sampleIndex][6]->Fill(cent, fPtProfilePi->GetBinCenter(i + 1), ((fPtProfilePi->GetBinContent(i + 1) / nSumEtaLeftPi) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][7]->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][8]->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        subSample[sampleIndex][9]->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
      }
    }

    if (nSumEtaRightHad > 0 && nSumEtaLeftHad > 0 && nSumEtaLeftKa > 0) {
      for (int i = 0; i < nbinsPid; i++) {
        histos.get<TProfile2D>(HIST("Prof_A_ka"))->Fill(cent, fPtProfileKa->GetBinCenter(i + 1), (fPtProfileKa->GetBinContent(i + 1) / nSumEtaLeftKa));
        histos.get<TProfile2D>(HIST("Prof_C_ka"))->Fill(cent, fPtProfileKa->GetBinCenter(i + 1), ((fPtProfileKa->GetBinContent(i + 1) / nSumEtaLeftKa) * (pTsumEtaRightHad / nSumEtaRightHad)));
        histos.get<TProfile2D>(HIST("Prof_Bone_ka"))->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        histos.get<TProfile2D>(HIST("Prof_Btwo_ka"))->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
        histos.get<TProfile2D>(HIST("Prof_D_ka"))->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));

        subSample[sampleIndex][10]->Fill(cent, fPtProfileKa->GetBinCenter(i + 1), (fPtProfileKa->GetBinContent(i + 1) / nSumEtaLeftKa));
        subSample[sampleIndex][11]->Fill(cent, fPtProfileKa->GetBinCenter(i + 1), ((fPtProfileKa->GetBinContent(i + 1) / nSumEtaLeftKa) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][12]->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][13]->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        subSample[sampleIndex][14]->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
      }
    }
    if (nSumEtaRightHad > 0 && nSumEtaLeftHad > 0 && nSumEtaLeftProt > 0) {
      for (int i = 1; i < nbinsPid; i++) {
        histos.get<TProfile2D>(HIST("Prof_A_prot"))->Fill(cent, fPtProfileProt->GetBinCenter(i + 1), (fPtProfileProt->GetBinContent(i + 1) / nSumEtaLeftProt));
        histos.get<TProfile2D>(HIST("Prof_C_prot"))->Fill(cent, fPtProfileProt->GetBinCenter(i + 1), ((fPtProfileProt->GetBinContent(i + 1) / nSumEtaLeftProt) * (pTsumEtaRightHad / nSumEtaRightHad)));
        histos.get<TProfile2D>(HIST("Prof_Bone_prot"))->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        histos.get<TProfile2D>(HIST("Prof_Btwo_prot"))->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
        histos.get<TProfile2D>(HIST("Prof_D_prot"))->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));

        subSample[sampleIndex][15]->Fill(cent, fPtProfileProt->GetBinCenter(i + 1), (fPtProfileProt->GetBinContent(i + 1) / nSumEtaLeftProt));
        subSample[sampleIndex][16]->Fill(cent, fPtProfileProt->GetBinCenter(i + 1), ((fPtProfileProt->GetBinContent(i + 1) / nSumEtaLeftProt) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][17]->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][18]->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        subSample[sampleIndex][19]->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
      }
    }

    if (nSumInWinA > cfgMinNoOfParticles && nSumInWinB > cfgMinNoOfParticles && nSumInWinC > cfgMinNoOfParticles) {
      double twoParCorr = (vecQInWinA * TComplex::Conjugate(vecQInWinC)).Re();
      twoParCorr *= 1.0 / (nSumInWinA * nSumInWinC);

      histos.get<TProfile2D>(HIST("Prof_XY"))->Fill(cent, 0.5, twoParCorr);
      subSampleV02[sampleIndex][0]->Fill(cent, 0.5, twoParCorr);

      if (cfgV02WeightedFill) {
        histos.get<TProfile2D>(HIST("Prof_XY_weighted"))->Fill(cent, 0.5, twoParCorr, (nSumInWinA * nSumInWinC));
        subSampleV02_weighted[sampleIndex][0]->Fill(cent, 0.5, twoParCorr, (nSumInWinA * nSumInWinC));
      }
      // hadrons
      for (int i = 0; i < cfgNbinsV02pt; i++) {
        double threeParCorrHad = (vecQInWinA * TComplex::Conjugate(vecQInWinC) * fPtProfileHadInWinB->GetBinContent(i + 1)).Re();
        threeParCorrHad *= 1.0 / (nSumInWinA * nSumInWinC * nSumInWinB);

        histos.get<TProfile2D>(HIST("Prof_XYZ_had"))->Fill(cent, fPtProfileHadInWinB->GetBinCenter(i + 1), threeParCorrHad);
        histos.get<TProfile2D>(HIST("Prof_Z_had"))->Fill(cent, fPtProfileHadInWinB->GetBinCenter(i + 1), (fPtProfileHadInWinB->GetBinContent(i + 1) / nSumInWinB));
        subSampleV02[sampleIndex][1]->Fill(cent, fPtProfileHadInWinB->GetBinCenter(i + 1), threeParCorrHad);
        subSampleV02[sampleIndex][2]->Fill(cent, fPtProfileHadInWinB->GetBinCenter(i + 1), (fPtProfileHadInWinB->GetBinContent(i + 1) / nSumInWinB));

        if (cfgV02WeightedFill) {
          histos.get<TProfile2D>(HIST("Prof_XYZ_weighted_had"))->Fill(cent, fPtProfileHadInWinB->GetBinCenter(i + 1), threeParCorrHad, (nSumInWinA * nSumInWinC));
          histos.get<TProfile2D>(HIST("Prof_Z_weighted_had"))->Fill(cent, fPtProfileHadInWinB->GetBinCenter(i + 1), (fPtProfileHadInWinB->GetBinContent(i + 1) / nSumInWinB));
          subSampleV02_weighted[sampleIndex][1]->Fill(cent, fPtProfileHadInWinB->GetBinCenter(i + 1), threeParCorrHad, (nSumInWinA * nSumInWinC));
          subSampleV02_weighted[sampleIndex][2]->Fill(cent, fPtProfileHadInWinB->GetBinCenter(i + 1), (fPtProfileHadInWinB->GetBinContent(i + 1) / nSumInWinB));
        }
      }

      // pions
      for (int i = 0; i < cfgNbinsV02pt; i++) {
        double threeParCorrPi = (vecQInWinA * TComplex::Conjugate(vecQInWinC) * fPtProfilePiInWinB->GetBinContent(i + 1)).Re();
        threeParCorrPi *= 1.0 / (nSumInWinA * nSumInWinC * nSumInWinB);

        histos.get<TProfile2D>(HIST("Prof_XYZ_pi"))->Fill(cent, fPtProfilePiInWinB->GetBinCenter(i + 1), threeParCorrPi);
        histos.get<TProfile2D>(HIST("Prof_Z_pi"))->Fill(cent, fPtProfilePiInWinB->GetBinCenter(i + 1), (fPtProfilePiInWinB->GetBinContent(i + 1) / nSumInWinB));
        subSampleV02[sampleIndex][3]->Fill(cent, fPtProfilePiInWinB->GetBinCenter(i + 1), threeParCorrPi);
        subSampleV02[sampleIndex][4]->Fill(cent, fPtProfilePiInWinB->GetBinCenter(i + 1), (fPtProfilePiInWinB->GetBinContent(i + 1) / nSumInWinB));

        if (cfgV02WeightedFill) {
          histos.get<TProfile2D>(HIST("Prof_XYZ_weighted_pi"))->Fill(cent, fPtProfilePiInWinB->GetBinCenter(i + 1), threeParCorrPi, (nSumInWinA * nSumInWinC));
          histos.get<TProfile2D>(HIST("Prof_Z_weighted_pi"))->Fill(cent, fPtProfilePiInWinB->GetBinCenter(i + 1), (fPtProfilePiInWinB->GetBinContent(i + 1) / nSumInWinB));
          subSampleV02_weighted[sampleIndex][3]->Fill(cent, fPtProfilePiInWinB->GetBinCenter(i + 1), threeParCorrPi, (nSumInWinA * nSumInWinC));
          subSampleV02_weighted[sampleIndex][4]->Fill(cent, fPtProfilePiInWinB->GetBinCenter(i + 1), (fPtProfilePiInWinB->GetBinContent(i + 1) / nSumInWinB));
        }
      }

      // kaons
      for (int i = 0; i < cfgNbinsV02pt; i++) {
        double threeParCorrKa = (vecQInWinA * TComplex::Conjugate(vecQInWinC) * fPtProfileKaInWinB->GetBinContent(i + 1)).Re();
        threeParCorrKa *= 1.0 / (nSumInWinA * nSumInWinC * nSumInWinB);
        histos.get<TProfile2D>(HIST("Prof_XYZ_ka"))->Fill(cent, fPtProfileKaInWinB->GetBinCenter(i + 1), threeParCorrKa);
        histos.get<TProfile2D>(HIST("Prof_Z_ka"))->Fill(cent, fPtProfileKaInWinB->GetBinCenter(i + 1), (fPtProfileKaInWinB->GetBinContent(i + 1) / nSumInWinB));

        subSampleV02[sampleIndex][5]->Fill(cent, fPtProfileKaInWinB->GetBinCenter(i + 1), threeParCorrKa);
        subSampleV02[sampleIndex][6]->Fill(cent, fPtProfileKaInWinB->GetBinCenter(i + 1), (fPtProfileKaInWinB->GetBinContent(i + 1) / nSumInWinB));

        if (cfgV02WeightedFill) {
          histos.get<TProfile2D>(HIST("Prof_XYZ_weighted_ka"))->Fill(cent, fPtProfileKaInWinB->GetBinCenter(i + 1), threeParCorrKa, (nSumInWinA * nSumInWinC));
          histos.get<TProfile2D>(HIST("Prof_Z_weighted_ka"))->Fill(cent, fPtProfileKaInWinB->GetBinCenter(i + 1), (fPtProfileKaInWinB->GetBinContent(i + 1) / nSumInWinB));
          subSampleV02_weighted[sampleIndex][3]->Fill(cent, fPtProfileKaInWinB->GetBinCenter(i + 1), threeParCorrKa, (nSumInWinA * nSumInWinC));
          subSampleV02_weighted[sampleIndex][4]->Fill(cent, fPtProfileKaInWinB->GetBinCenter(i + 1), (fPtProfileKaInWinB->GetBinContent(i + 1) / nSumInWinB));
        }
      }

      // protons
      for (int i = 1; i < cfgNbinsV02pt; i++) {
        double threeParCorrProt = (vecQInWinA * TComplex::Conjugate(vecQInWinC) * fPtProfileProtInWinB->GetBinContent(i + 1)).Re();
        threeParCorrProt *= 1.0 / (nSumInWinA * nSumInWinC * nSumInWinB);
        histos.get<TProfile2D>(HIST("Prof_XYZ_prot"))->Fill(cent, fPtProfileProtInWinB->GetBinCenter(i + 1), threeParCorrProt);
        histos.get<TProfile2D>(HIST("Prof_Z_prot"))->Fill(cent, fPtProfileProtInWinB->GetBinCenter(i + 1), (fPtProfileProtInWinB->GetBinContent(i + 1) / nSumInWinB));

        subSampleV02[sampleIndex][7]->Fill(cent, fPtProfileProtInWinB->GetBinCenter(i + 1), threeParCorrProt);
        subSampleV02[sampleIndex][8]->Fill(cent, fPtProfileProtInWinB->GetBinCenter(i + 1), (fPtProfileProtInWinB->GetBinContent(i + 1) / nSumInWinB));

        if (cfgV02WeightedFill) {
          histos.get<TProfile2D>(HIST("Prof_XYZ_weighted_prot"))->Fill(cent, fPtProfileProtInWinB->GetBinCenter(i + 1), threeParCorrProt, (nSumInWinA * nSumInWinC));
          histos.get<TProfile2D>(HIST("Prof_Z_weighted_prot"))->Fill(cent, fPtProfileProtInWinB->GetBinCenter(i + 1), (fPtProfileProtInWinB->GetBinContent(i + 1) / nSumInWinB));
          subSampleV02_weighted[sampleIndex][3]->Fill(cent, fPtProfileProtInWinB->GetBinCenter(i + 1), threeParCorrProt, (nSumInWinA * nSumInWinC));
          subSampleV02_weighted[sampleIndex][4]->Fill(cent, fPtProfileProtInWinB->GetBinCenter(i + 1), (fPtProfileProtInWinB->GetBinContent(i + 1) / nSumInWinB));
        }
      }
    }

    fPtProfileHad->Delete();
    fPtProfilePi->Delete();
    fPtProfileKa->Delete();
    fPtProfileProt->Delete();

    fPtProfileHadInWinB->Delete();
    fPtProfilePiInWinB->Delete();
    fPtProfileKaInWinB->Delete();
    fPtProfileProtInWinB->Delete();

  } // End process loop
  PROCESS_SWITCH(V0ptHadPiKaProt, processData, "Process Real Data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<V0ptHadPiKaProt>(cfgc)};
  return workflow;
}
