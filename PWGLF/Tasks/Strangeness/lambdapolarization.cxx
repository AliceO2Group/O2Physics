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

/// \author Junlee Kim (jikim1290@gmail.com)

#include <cmath>
#include <array>
#include <cstdlib>
#include <chrono>
#include <string>
#include <vector>

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TVector2.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include <TMath.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/PIDResponseITS.h"

#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"

#include "CommonConstants/PhysicsConstants.h"

#include "ReconstructionDataFormats/Track.h"

#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGMM/Mult/DataModel/Index.h" // for Particles2Tracks table

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct lambdapolarization {
  //  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults, aod::Qvectors>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using V0TrackCandidate = aod::V0Datas;

  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL",
                                     "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than",
                                      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(),
                                      "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  Configurable<float> cfgCentSel{"cfgCentSel", 80., "Centrality selection"};
  Configurable<int> cfgCentEst{"cfgCentEst", 1, "Centrality estimator, 1: FT0C, 2: FT0M"};

  Configurable<bool> cfgPVSel{"cfgPVSel", false, "Additional PV selection flag for syst"};
  Configurable<float> cfgPV{"cfgPV", 8.0, "Additional PV selection range for syst"};
  Configurable<bool> cfgAddEvtSelPileup{"cfgAddEvtSelPileup", false, "flag for additional pileup selection"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "maximum occupancy of tracks in neighbouring collisions in a given time range"};

  Configurable<float> cfgv0radiusMin{"cfgv0radiusMin", 1.2, "minimum decay radius"};
  Configurable<float> cfgDCAPrToPVMin{"cfgDCAPrToPVMin", 0.05, "minimum DCA to PV for proton track"};
  Configurable<float> cfgDCAPiToPVMin{"cfgDCAPiToPVMin", 0.1, "minimum DCA to PV for pion track"};
  Configurable<float> cfgv0CosPA{"cfgv0CosPA", 0.995, "minimum v0 cosine"};
  Configurable<float> cfgDCAV0Dau{"cfgDCAV0Dau", 1.0, "maximum DCA between daughters"};

  Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0, "minimum pT for lambda"};
  Configurable<float> cfgV0EtaMin{"cfgV0EtaMin", -0.5, "maximum rapidity"};
  Configurable<float> cfgV0EtaMax{"cfgV0EtaMax", 0.5, "maximum rapidity"};
  Configurable<float> cfgV0LifeTime{"cfgV0LifeTime", 30., "maximum lambda lifetime"};

  Configurable<bool> cfgQAv0{"cfgQAv0", false, "QA plot"};

  Configurable<int> cfgDaughTPCnclsMin{"cfgDaughTPCnclsMin", 70, "minimum fired crossed rows"};
  Configurable<float> cfgDaughPIDCutsTPCPr{"cfgDaughPIDCutsTPCPr", 5, "proton nsigma for TPC"};
  Configurable<float> cfgDaughPIDCutsTPCPi{"cfgDaughPIDCutsTPCPi", 5, "pion nsigma for TPC"};
  Configurable<float> cfgDaughEtaMin{"cfgDaughEtaMin", -0.8, "minimum daughter eta"};
  Configurable<float> cfgDaughEtaMax{"cfgDaughEtaMax", 0.8, "maximum daughter eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.5, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.5, "minimum daughter pion pt"};

  Configurable<int> cfgnMods{"cfgnMods", 1, "The number of modulations of interest starting from 2"};
  Configurable<int> cfgNQvec{"cfgNQvec", 7, "The number of total Qvectors for looping over the task"};

  Configurable<std::string> cfgQvecDetName{"cfgQvecDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgQvecRefAName{"cfgQvecRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgQvecRefBName{"cfgQvecRefBName", "TPCneg", "The name of detector for reference B"};

  Configurable<bool> cfgPhiDepStudy{"cfgPhiDepStudy", false, "cfg for phi dependent study"};
  Configurable<bool> cfgUSESP{"cfgUSESP", false, "cfg for sp"};
  Configurable<float> cfgPhiDepSig{"cfgPhiDepSig", 0.2, "cfg for significance on phi dependent study"};

  Configurable<bool> cfgShiftCorr{"cfgShiftCorr", false, "additional shift correction"};
  Configurable<bool> cfgShiftCorrDef{"cfgShiftCorrDef", false, "additional shift correction definition"};
  Configurable<std::string> cfgShiftPath{"cfgShiftPath", "Users/j/junlee/Qvector/QvecCalib/Shift", "Path for Shift"};

  Configurable<bool> cfgEffCor{"cfgEffCor", false, "flag to apply efficiency correction"};
  Configurable<std::string> cfgEffCorPath{"cfgEffCorPath", "", "path for pseudo efficiency correction"};

  Configurable<bool> cfgAccCor{"cfgAccCor", false, "flag to apply acceptance correction"};
  Configurable<std::string> cfgAccCorPath{"cfgAccCorPath", "", "path for pseudo acceptance correction"};

  Configurable<bool> cfgCalcCum{"cfgCalcCum", false, "flag to calculate cumulants of cossin"};
  Configurable<bool> cfgCalcCum1{"cfgCalcCum1", false, "flag to calculate cumulants of coscos"};

  Configurable<bool> cfgRapidityDep{"cfgRapidityDep", false, "flag for rapidity dependent study"};
  Configurable<bool> cfgAccAzimuth{"cfgAccAzimuth", false, "flag for azimuth closure study"};

  ConfigurableAxis massAxis{"massAxis", {30, 1.1, 1.13}, "Invariant mass axis"};
  ConfigurableAxis ptAxis{"ptAxis", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "Transverse momentum bins"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100}, "Centrality interval"};
  ConfigurableAxis cosAxis{"cosAxis", {110, -1.05, 1.05}, "Cosine axis"};
  ConfigurableAxis RapAxis{"RapAxis", {10, -0.5, 0.5}, "Rapidity axis"};
  ConfigurableAxis qqAxis{"qqAxis", {100, -0.1, 0.1}, "qq axis"};

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;

  int DetId;
  int RefAId;
  int RefBId;

  int QvecDetInd;
  int QvecRefAInd;
  int QvecRefBInd;

  float centrality;

  double angle;
  double psi;
  double relphi;

  int currentRunNumber = -999;
  int lastRunNumber = -999;
  std::vector<TProfile3D*> shiftprofile{};
  TProfile2D* EffMap = nullptr;
  TProfile2D* AccMap = nullptr;

  std::string fullCCDBShiftCorrPath;

  double GetPhiInRange(double phi)
  {
    double result = phi;
    while (result < 0) {
      result = result + 2. * TMath::Pi() / 2;
    }
    while (result > 2. * TMath::Pi() / 2) {
      result = result - 2. * TMath::Pi() / 2;
    }
    return result;
  }

  template <typename T>
  int GetDetId(const T& name)
  {
    if (name.value == "FT0C") {
      return 0;
    } else if (name.value == "FT0A") {
      return 1;
    } else if (name.value == "FT0M") {
      return 2;
    } else if (name.value == "FV0A") {
      return 3;
    } else if (name.value == "TPCpos") {
      return 4;
    } else if (name.value == "TPCneg") {
      return 5;
    } else if (name.value == "TPCall") {
      return 6;
    } else {
      return 0;
    }
  }

  void init(o2::framework::InitContext&)
  {
    AxisSpec centQaAxis = {80, 0.0, 80.0};
    AxisSpec PVzQaAxis = {300, -15.0, 15.0};
    AxisSpec epAxis = {6, 0.0, 2.0 * constants::math::PI};
    AxisSpec epQaAxis = {100, -1.0 * constants::math::PI, constants::math::PI};

    AxisSpec pidAxis = {100, -10, 10};

    AxisSpec shiftAxis = {10, 0, 10, "shift"};
    AxisSpec basisAxis = {20, 0, 20, "basis"};

    for (auto i = 2; i < cfgnMods + 2; i++) {
      histos.add(Form("psi%d/h_lambda_cos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      histos.add(Form("psi%d/h_alambda_cos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      histos.add(Form("psi%d/h_lambda_cos2", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      histos.add(Form("psi%d/h_alambda_cos2", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});

      if (cfgRapidityDep) {
        histos.add(Form("psi%d/h_lambda_cos2_rap", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, RapAxis}});
        histos.add(Form("psi%d/h_alambda_cos2_rap", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, RapAxis}});
      }

      histos.add(Form("psi%d/h_lambda_cossin", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add(Form("psi%d/h_alambda_cossin", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      if (cfgAccAzimuth) {
        histos.add(Form("psi%d/h_lambda_coscos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
        histos.add(Form("psi%d/h_alambda_coscos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      }

      histos.add(Form("psi%d/h_lambda_vncos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add(Form("psi%d/h_lambda_vnsin", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add(Form("psi%d/h_alambda_vncos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add(Form("psi%d/h_alambda_vnsin", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
    }
    histos.add("QA/ptspec_l", "", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});
    histos.add("QA/ptspec_al", "", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});
    histos.add("QA/ptspecCor_l", "", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});
    histos.add("QA/ptspecCor_al", "", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});

    if (cfgCalcCum) {
      histos.add("psi2/QA/cosTheta_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/sinPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/sinPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosTheta_cosPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_cosPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosTheta_sinPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_sinPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosPhi_sinPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/sinPhi_cosPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosTheta_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/sinPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/sinPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosTheta_cosPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_cosPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosTheta_sinPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_sinPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosPhi_sinPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/sinPhi_cosPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
    }

    if (cfgCalcCum1) {
      histos.add("psi2/QA/cosTheta_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosPhi_cosPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_cosPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_cosPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/sinPhi_sinPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_sinPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_sinPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/sinPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/sinPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosTheta_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosPhi_cosPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_cosPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_cosPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/sinPhi_sinPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_sinPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_sinPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/sinPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/sinPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
    }

    if (cfgQAv0) {
      histos.add("QA/CentDist", "", {HistType::kTH1F, {centQaAxis}});
      histos.add("QA/PVzDist", "", {HistType::kTH1F, {PVzQaAxis}});

      histos.add("QA/nsigma_tpc_pt_ppr", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_ppi", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_mpr", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_mpi", "", {HistType::kTH2F, {ptAxis, pidAxis}});

      for (auto i = 2; i < cfgnMods + 2; i++) {
        histos.add(Form("psi%d/QA/EP_Det", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
        histos.add(Form("psi%d/QA/EP_RefA", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
        histos.add(Form("psi%d/QA/EP_RefB", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});

        histos.add(Form("psi%d/QA/qqAxis_Det_RefA_xx", i), "", {HistType::kTH2F, {centQaAxis, qqAxis}});
        histos.add(Form("psi%d/QA/qqAxis_Det_RefB_xx", i), "", {HistType::kTH2F, {centQaAxis, qqAxis}});
        histos.add(Form("psi%d/QA/qqAxis_RefA_RefB_xx", i), "", {HistType::kTH2F, {centQaAxis, qqAxis}});

        histos.add(Form("psi%d/QA/qqAxis_Det_RefA_yy", i), "", {HistType::kTH2F, {centQaAxis, qqAxis}});
        histos.add(Form("psi%d/QA/qqAxis_Det_RefB_yy", i), "", {HistType::kTH2F, {centQaAxis, qqAxis}});
        histos.add(Form("psi%d/QA/qqAxis_RefA_RefB_yy", i), "", {HistType::kTH2F, {centQaAxis, qqAxis}});

        histos.add(Form("psi%d/QA/EPRes_Det_RefA", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EPRes_Det_RefB", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EPRes_RefA_RefB", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});

        histos.add(Form("psi%d/QA/EP_FT0C_shifted", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
        histos.add(Form("psi%d/QA/EP_FT0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
        histos.add(Form("psi%d/QA/EP_FV0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});

        histos.add(Form("psi%d/QA/EPRes_FT0C_FT0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EPRes_FT0C_FV0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EPRes_FT0A_FV0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
      }
    }

    if (doprocessMC_ITSTPC) {
      histos.add("hImpactParameter", "Impact parameter", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("hEventPlaneAngle", "hEventPlaneAngle", kTH1F, {{200, -2.0f * TMath::Pi(), 2.0f * TMath::Pi()}});
      histos.add("hEventPlaneAngleRec", "hEventPlaneAngleRec", kTH1F, {{200, -2.0f * TMath::Pi(), 2.0f * TMath::Pi()}});
      histos.add("hNchVsImpactParameter", "hNchVsImpactParameter", kTH2F, {{200, 0.0f, 20.0f}, {500, -0.5f, 5000.5f}});
      histos.add("hSparseMCGenWeight", "hSparseMCGenWeight", HistType::kTHnSparseF, {centAxis, {36, 0.0f, TMath::Pi()}, {50, 0.0f, 1}, ptAxis, {8, -0.8, 0.8}});
      histos.add("hSparseMCRecWeight", "hSparseMCRecWeight", HistType::kTHnSparseF, {centAxis, {36, 0.0f, TMath::Pi()}, {50, 0.0f, 1}, ptAxis, {8, -0.8, 0.8}});
      histos.add("hSparseMCRecAllTrackWeight", "hSparseMCRecAllTrackWeight", HistType::kTHnSparseF, {centAxis, {36, 0.0, TMath::Pi()}, {50, 0.0f, 1}, ptAxis, {8, -0.8, 0.8}});
    }

    if (cfgShiftCorrDef) {
      for (auto i = 2; i < cfgnMods + 2; i++) {
        histos.add(Form("psi%d/ShiftFIT", i), "", kTProfile3D, {centQaAxis, basisAxis, shiftAxis});
      }
    }

    DetId = GetDetId(cfgQvecDetName);
    RefAId = GetDetId(cfgQvecRefAName);
    RefBId = GetDetId(cfgQvecRefBName);

    if (DetId == RefAId || DetId == RefBId || RefAId == RefBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      DetId = 0;
      RefAId = 4;
      RefBId = 5;
    }

    fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
    fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);

    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  ROOT::Math::PxPyPzMVector ProtonVec, PionVec, LambdaVec, ProtonBoostedVec, PionBoostedVec;

  template <typename TCollision>
  bool eventSelected(TCollision collision)
  {
    if (!collision.sel8()) {
      return 0;
    }

    if (cfgCentSel < centrality) {
      return 0;
    }
    /*
        auto multNTracksPV = collision.multNTracksPV();
        if (multNTracksPV < fMultPVCutLow->Eval(centrality)) {
          return 0;
        }
        if (multNTracksPV > fMultPVCutHigh->Eval(centrality)) {
          return 0;
        }
    */
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    if (cfgPVSel && std::abs(collision.posZ()) > cfgPV) {
      return 0;
    }
    if (cfgAddEvtSelPileup && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
    }
    if (collision.trackOccupancyInTimeRange() > cfgMaxOccupancy || collision.trackOccupancyInTimeRange() < cfgMinOccupancy) {
      return 0;
    }

    return 1;
  } // event selection

  template <typename TCollision, typename V0>
  bool SelectionV0(TCollision const& collision, V0 const& candidate, int LambdaTag)
  {
    if (candidate.v0radius() < cfgv0radiusMin)
      return false;
    if (LambdaTag) {
      if (std::abs(candidate.dcapostopv()) < cfgDCAPrToPVMin)
        return false;
      if (std::abs(candidate.dcanegtopv()) < cfgDCAPiToPVMin)
        return false;
    } else if (!LambdaTag) {
      if (std::abs(candidate.dcapostopv()) < cfgDCAPiToPVMin)
        return false;
      if (std::abs(candidate.dcanegtopv()) < cfgDCAPrToPVMin)
        return false;
    }
    if (candidate.v0cosPA() < cfgv0CosPA)
      return false;
    if (std::abs(candidate.dcaV0daughters()) > cfgDCAV0Dau)
      return false;
    if (candidate.pt() < cfgV0PtMin)
      return false;
    if (candidate.yLambda() < cfgV0EtaMin)
      return false;
    if (candidate.yLambda() > cfgV0EtaMax)
      return false;
    if (candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda > cfgV0LifeTime)
      return false;

    return true;
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, int pid) // pid 0: proton, pid 1: pion
  {
    if (track.tpcNClsFound() < cfgDaughTPCnclsMin)
      return false;
    if (pid == 0 && std::abs(track.tpcNSigmaPr()) > cfgDaughPIDCutsTPCPr)
      return false;
    if (pid == 1 && std::abs(track.tpcNSigmaPi()) > cfgDaughPIDCutsTPCPi)
      return false;
    if (track.eta() > cfgDaughEtaMax)
      return false;
    if (track.eta() < cfgDaughEtaMin)
      return false;
    if (pid == 0 && track.pt() < cfgDaughPrPt)
      return false;
    if (pid == 1 && track.pt() < cfgDaughPiPt)
      return false;

    return true;
  }

  template <typename TCollision>
  void FillShiftCorrection(TCollision const& collision, int nmode)
  {
    QvecDetInd = DetId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefAInd = RefAId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefBInd = RefBId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    for (int ishift = 1; ishift <= 10; ishift++) {
      if (nmode == 2) {
        histos.fill(HIST("psi2/ShiftFIT"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi2/ShiftFIT"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi2/ShiftFIT"), centrality, 2.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi2/ShiftFIT"), centrality, 3.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi2/ShiftFIT"), centrality, 4.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi2/ShiftFIT"), centrality, 5.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode)));
      } else if (nmode == 3) {
        histos.fill(HIST("psi3/ShiftFIT"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi3/ShiftFIT"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi3/ShiftFIT"), centrality, 2.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi3/ShiftFIT"), centrality, 3.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi3/ShiftFIT"), centrality, 4.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi3/ShiftFIT"), centrality, 5.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode)));
      } else if (nmode == 4) {
        histos.fill(HIST("psi4/ShiftFIT"), centrality, 0.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi4/ShiftFIT"), centrality, 1.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi4/ShiftFIT"), centrality, 2.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi4/ShiftFIT"), centrality, 3.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi4/ShiftFIT"), centrality, 4.5, ishift - 0.5, TMath::Sin(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi4/ShiftFIT"), centrality, 5.5, ishift - 0.5, TMath::Cos(ishift * static_cast<float>(nmode) * TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode)));
      }
    }
  }

  template <typename TCollision>
  void FillEPQA(TCollision const& collision, int nmode)
  {
    QvecDetInd = DetId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefAInd = RefAId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefBInd = RefBId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;

    if (collision.qvecAmp()[DetId] < 1e-5 || collision.qvecAmp()[RefAId] < 1e-5 || collision.qvecAmp()[RefBId] < 1e-5)
      return;

    if (nmode == 2) {
      histos.fill(HIST("psi2/QA/EP_Det"), centrality, TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi2/QA/EP_RefA"), centrality, TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi2/QA/EP_RefB"), centrality, TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode));

      histos.fill(HIST("psi2/QA/qqAxis_Det_RefA_xx"), centrality, collision.qvecRe()[QvecDetInd] * collision.qvecRe()[QvecRefAInd]);
      histos.fill(HIST("psi2/QA/qqAxis_Det_RefB_xx"), centrality, collision.qvecRe()[QvecDetInd] * collision.qvecRe()[QvecRefBInd]);
      histos.fill(HIST("psi2/QA/qqAxis_RefA_RefB_xx"), centrality, collision.qvecRe()[QvecRefAInd] * collision.qvecRe()[QvecRefBInd]);

      histos.fill(HIST("psi2/QA/qqAxis_Det_RefA_yy"), centrality, collision.qvecIm()[QvecDetInd] * collision.qvecIm()[QvecRefAInd]);
      histos.fill(HIST("psi2/QA/qqAxis_Det_RefB_yy"), centrality, collision.qvecIm()[QvecDetInd] * collision.qvecIm()[QvecRefBInd]);
      histos.fill(HIST("psi2/QA/qqAxis_RefA_RefB_yy"), centrality, collision.qvecIm()[QvecRefAInd] * collision.qvecIm()[QvecRefBInd]);

      histos.fill(HIST("psi2/QA/EPRes_Det_RefA"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) - TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd])));
      histos.fill(HIST("psi2/QA/EPRes_Det_RefB"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) - TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd])));
      histos.fill(HIST("psi2/QA/EPRes_RefA_RefB"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) - TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd])));
    } else if (nmode == 3) {
      histos.fill(HIST("psi3/QA/EP_Det"), centrality, TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi3/QA/EP_RefA"), centrality, TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi3/QA/EP_RefB"), centrality, TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode));

      histos.fill(HIST("psi3/QA/qqAxis_Det_RefA_xx"), centrality, collision.qvecRe()[QvecDetInd] * collision.qvecRe()[QvecRefAInd]);
      histos.fill(HIST("psi3/QA/qqAxis_Det_RefB_xx"), centrality, collision.qvecRe()[QvecDetInd] * collision.qvecRe()[QvecRefBInd]);
      histos.fill(HIST("psi3/QA/qqAxis_RefA_RefB_xx"), centrality, collision.qvecRe()[QvecRefAInd] * collision.qvecRe()[QvecRefBInd]);

      histos.fill(HIST("psi3/QA/qqAxis_Det_RefA_yy"), centrality, collision.qvecIm()[QvecDetInd] * collision.qvecIm()[QvecRefAInd]);
      histos.fill(HIST("psi3/QA/qqAxis_Det_RefB_yy"), centrality, collision.qvecIm()[QvecDetInd] * collision.qvecIm()[QvecRefBInd]);
      histos.fill(HIST("psi3/QA/qqAxis_RefA_RefB_yy"), centrality, collision.qvecIm()[QvecRefAInd] * collision.qvecIm()[QvecRefBInd]);

      histos.fill(HIST("psi3/QA/EPRes_Det_RefA"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) - TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd])));
      histos.fill(HIST("psi3/QA/EPRes_Det_RefB"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) - TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd])));
      histos.fill(HIST("psi3/QA/EPRes_RefA_RefB"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) - TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd])));
    } else if (nmode == 4) {
      histos.fill(HIST("psi4/QA/EP_Det"), centrality, TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi4/QA/EP_RefA"), centrality, TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi4/QA/EP_RefB"), centrality, TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode));

      histos.fill(HIST("psi4/QA/qqAxis_Det_RefA_xx"), centrality, collision.qvecRe()[QvecDetInd] * collision.qvecRe()[QvecRefAInd]);
      histos.fill(HIST("psi4/QA/qqAxis_Det_RefB_xx"), centrality, collision.qvecRe()[QvecDetInd] * collision.qvecRe()[QvecRefBInd]);
      histos.fill(HIST("psi4/QA/qqAxis_RefA_RefB_xx"), centrality, collision.qvecRe()[QvecRefAInd] * collision.qvecRe()[QvecRefBInd]);

      histos.fill(HIST("psi4/QA/qqAxis_Det_RefA_yy"), centrality, collision.qvecIm()[QvecDetInd] * collision.qvecIm()[QvecRefAInd]);
      histos.fill(HIST("psi4/QA/qqAxis_Det_RefB_yy"), centrality, collision.qvecIm()[QvecDetInd] * collision.qvecIm()[QvecRefBInd]);
      histos.fill(HIST("psi4/QA/qqAxis_RefA_RefB_yy"), centrality, collision.qvecIm()[QvecRefAInd] * collision.qvecIm()[QvecRefBInd]);

      histos.fill(HIST("psi4/QA/EPRes_Det_RefA"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) - TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd])));
      histos.fill(HIST("psi4/QA/EPRes_Det_RefB"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) - TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd])));
      histos.fill(HIST("psi4/QA/EPRes_RefA_RefB"), centrality, TMath::Cos(TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) - TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd])));
    }

    if (cfgShiftCorr) {
      auto deltapsiFT0C = 0.0;
      auto deltapsiFT0A = 0.0;
      auto deltapsiFV0A = 0.0;

      auto psidefFT0C = TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode);
      auto psidefFT0A = TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode);
      auto psidefFV0A = TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode);
      for (int ishift = 1; ishift <= 10; ishift++) {
        auto coeffshiftxFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 0.5, ishift - 0.5));
        auto coeffshiftyFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 1.5, ishift - 0.5));
        auto coeffshiftxFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 2.5, ishift - 0.5));
        auto coeffshiftyFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 3.5, ishift - 0.5));
        auto coeffshiftxFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 4.5, ishift - 0.5));
        auto coeffshiftyFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 5.5, ishift - 0.5));

        deltapsiFT0C += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0C * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFT0C) + coeffshiftyFT0C * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFT0C)));
        deltapsiFT0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0A * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFT0A) + coeffshiftyFT0A * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFT0A)));
        deltapsiFV0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFV0A * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFV0A) + coeffshiftyFV0A * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFV0A)));
      }
      if (nmode == 2) {
        histos.fill(HIST("psi2/QA/EP_FT0C_shifted"), centrality, psidefFT0C + deltapsiFT0C);
        histos.fill(HIST("psi2/QA/EP_FT0A_shifted"), centrality, psidefFT0A + deltapsiFT0A);
        histos.fill(HIST("psi2/QA/EP_FV0A_shifted"), centrality, psidefFV0A + deltapsiFV0A);

        histos.fill(HIST("psi2/QA/EPRes_FT0C_FT0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFT0A - deltapsiFT0A)));
        histos.fill(HIST("psi2/QA/EPRes_FT0C_FV0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFV0A - deltapsiFV0A)));
        histos.fill(HIST("psi2/QA/EPRes_FT0A_FV0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0A + deltapsiFT0A - psidefFV0A - deltapsiFV0A)));
      } else if (nmode == 3) {
        histos.fill(HIST("psi3/QA/EP_FT0C_shifted"), centrality, psidefFT0C + deltapsiFT0C);
        histos.fill(HIST("psi3/QA/EP_FT0A_shifted"), centrality, psidefFT0A + deltapsiFT0A);
        histos.fill(HIST("psi3/QA/EP_FV0A_shifted"), centrality, psidefFV0A + deltapsiFV0A);

        histos.fill(HIST("psi3/QA/EPRes_FT0C_FT0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFT0A - deltapsiFT0A)));
        histos.fill(HIST("psi3/QA/EPRes_FT0C_FV0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFV0A - deltapsiFV0A)));
        histos.fill(HIST("psi3/QA/EPRes_FT0A_FV0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0A + deltapsiFT0A - psidefFV0A - deltapsiFV0A)));
      } else if (nmode == 4) {
        histos.fill(HIST("psi4/QA/EP_FT0C_shifted"), centrality, psidefFT0C + deltapsiFT0C);
        histos.fill(HIST("psi4/QA/EP_FT0A_shifted"), centrality, psidefFT0A + deltapsiFT0A);
        histos.fill(HIST("psi4/QA/EP_FV0A_shifted"), centrality, psidefFV0A + deltapsiFV0A);

        histos.fill(HIST("psi4/QA/EPRes_FT0C_FT0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFT0A - deltapsiFT0A)));
        histos.fill(HIST("psi4/QA/EPRes_FT0C_FV0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFV0A - deltapsiFV0A)));
        histos.fill(HIST("psi4/QA/EPRes_FT0A_FV0A_shifted"), centrality, TMath::Cos(static_cast<float>(nmode) * (psidefFT0A + deltapsiFT0A - psidefFV0A - deltapsiFV0A)));
      }
    }
  }

  template <typename TCollision, typename V0>
  void FillHistograms(TCollision const& collision, V0 const& V0s, int nmode)
  {
    QvecDetInd = DetId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefAInd = RefAId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefBInd = RefBId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;

    for (auto& v0 : V0s) {
      auto postrack = v0.template posTrack_as<TrackCandidates>();
      auto negtrack = v0.template negTrack_as<TrackCandidates>();

      double nTPCSigmaPosPr = postrack.tpcNSigmaPr();
      double nTPCSigmaNegPi = negtrack.tpcNSigmaPi();

      double nTPCSigmaNegPr = negtrack.tpcNSigmaPr();
      double nTPCSigmaPosPi = postrack.tpcNSigmaPi();

      if (cfgQAv0 && nmode == 2) {
        histos.fill(HIST("QA/nsigma_tpc_pt_ppr"), postrack.pt(), nTPCSigmaPosPr);
        histos.fill(HIST("QA/nsigma_tpc_pt_ppi"), postrack.pt(), nTPCSigmaPosPi);

        histos.fill(HIST("QA/nsigma_tpc_pt_mpr"), negtrack.pt(), nTPCSigmaNegPr);
        histos.fill(HIST("QA/nsigma_tpc_pt_mpi"), negtrack.pt(), nTPCSigmaNegPi);
      }

      int LambdaTag = 0;
      int aLambdaTag = 0;

      if (isSelectedV0Daughter(postrack, 0) && isSelectedV0Daughter(negtrack, 1)) {
        LambdaTag = 1;
      }
      if (isSelectedV0Daughter(negtrack, 0) && isSelectedV0Daughter(postrack, 1)) {
        aLambdaTag = 1;
      }

      if (LambdaTag == aLambdaTag)
        continue;

      if (!SelectionV0(collision, v0, LambdaTag))
        continue;

      if (LambdaTag) {
        ProtonVec = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPr);
        PionVec = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPi);
      }
      if (aLambdaTag) {
        ProtonVec = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPr);
        PionVec = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPi);
      }
      LambdaVec = ProtonVec + PionVec;
      LambdaVec.SetM(massLambda);

      ROOT::Math::Boost boost{LambdaVec.BoostToCM()};
      ProtonBoostedVec = boost(ProtonVec);

      angle = ProtonBoostedVec.Pz() / ProtonBoostedVec.P();
      psi = TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode);
      relphi = TVector2::Phi_0_2pi(static_cast<float>(nmode) * (LambdaVec.Phi() - psi));

      if (cfgShiftCorr) {
        auto deltapsiFT0C = 0.0;
        auto deltapsiFT0A = 0.0;
        auto deltapsiFV0A = 0.0;

        auto psidefFT0C = TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode);
        auto psidefFT0A = TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode);
        auto psidefFV0A = TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode);
        for (int ishift = 1; ishift <= 10; ishift++) {
          auto coeffshiftxFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 0.5, ishift - 0.5));
          auto coeffshiftyFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 1.5, ishift - 0.5));
          auto coeffshiftxFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 2.5, ishift - 0.5));
          auto coeffshiftyFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 3.5, ishift - 0.5));
          auto coeffshiftxFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 4.5, ishift - 0.5));
          auto coeffshiftyFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 5.5, ishift - 0.5));

          deltapsiFT0C += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0C * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFT0C) + coeffshiftyFT0C * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFT0C)));
          deltapsiFT0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0A * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFT0A) + coeffshiftyFT0A * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFT0A)));
          deltapsiFV0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFV0A * TMath::Cos(ishift * static_cast<float>(nmode) * psidefFV0A) + coeffshiftyFV0A * TMath::Sin(ishift * static_cast<float>(nmode) * psidefFV0A)));
        }
        psi += deltapsiFT0C;
        relphi = TVector2::Phi_0_2pi(static_cast<float>(nmode) * (LambdaVec.Phi() - psidefFT0C - deltapsiFT0C));
      }

      if (cfgPhiDepStudy && cfgPhiDepSig * std::abs(TMath::Sin(relphi)) > gRandom->Uniform(0, 1)) {
        continue;
      }

      if (LambdaTag) {
        histos.fill(HIST("QA/ptspec_l"), v0.mLambda(), v0.pt(), centrality);
        if (cfgEffCor) {
          histos.fill(HIST("QA/ptspecCor_l"), v0.mLambda(), v0.pt(), centrality,
                      1.0 / EffMap->GetBinContent(EffMap->GetXaxis()->FindBin(v0.pt()), EffMap->GetYaxis()->FindBin(centrality)));
        }
      }
      if (aLambdaTag) {
        histos.fill(HIST("QA/ptspec_al"), v0.mAntiLambda(), v0.pt(), centrality);
        if (cfgEffCor) {
          histos.fill(HIST("QA/ptspecCor_al"), v0.mAntiLambda(), v0.pt(), centrality,
                      1.0 / EffMap->GetBinContent(EffMap->GetXaxis()->FindBin(v0.pt()), EffMap->GetYaxis()->FindBin(centrality)));
        }
      }
      double weight = 1.0;
      weight *= cfgEffCor ? 1.0 / EffMap->GetBinContent(EffMap->GetXaxis()->FindBin(v0.pt()), EffMap->GetYaxis()->FindBin(centrality)) : 1.;
      weight *= cfgAccCor ? 1.0 / AccMap->GetBinContent(AccMap->GetXaxis()->FindBin(v0.pt()), AccMap->GetYaxis()->FindBin(v0.yLambda())) : 1.;

      double qvecMag = 1.0;
      if (cfgUSESP)
        qvecMag *= TMath::Sqrt(TMath::Power(collision.qvecIm()[3 + (nmode - 2) * 28], 2) + TMath::Power(collision.qvecRe()[3 + (nmode - 2) * 28], 2));

      if (nmode == 2) { ////////////
        if (LambdaTag) {
          histos.fill(HIST("psi2/h_lambda_cos"), v0.mLambda(), v0.pt(), angle * weight, centrality, relphi);
          histos.fill(HIST("psi2/h_lambda_cos2"), v0.mLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi2/h_lambda_cossin"), v0.mLambda(), v0.pt(), angle * TMath::Sin(relphi) * weight, centrality);
          histos.fill(HIST("psi2/h_lambda_vncos"), v0.mLambda(), v0.pt(), qvecMag * TMath::Cos(relphi) * weight, centrality);
          histos.fill(HIST("psi2/h_lambda_vnsin"), v0.mLambda(), v0.pt(), TMath::Sin(relphi), centrality);

          if (cfgRapidityDep) {
            histos.fill(HIST("psi2/h_lambda_cos2_rap"), v0.mLambda(), v0.pt(), angle * angle, centrality, v0.yLambda(), weight);
          }

          if (cfgAccAzimuth) {
            histos.fill(HIST("psi2/h_lambda_coscos"), v0.mLambda(), v0.pt(), angle * TMath::Cos(relphi), centrality, weight);
          }

          if (cfgCalcCum) {
            histos.fill(HIST("psi2/QA/cosTheta_l"), v0.mLambda(), v0.pt(), angle, centrality);
            histos.fill(HIST("psi2/QA/cosPsi_l"), v0.mLambda(), v0.pt(), TMath::Cos(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosPhi_l"), v0.mLambda(), v0.pt(), TMath::Cos(v0.phi() * 2.0), centrality);

            histos.fill(HIST("psi2/QA/sinPsi_l"), v0.mLambda(), v0.pt(), TMath::Sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/sinPhi_l"), v0.mLambda(), v0.pt(), TMath::Sin(v0.phi() * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosTheta_cosPhi_l"), v0.mLambda(), v0.pt(), angle * TMath::Cos(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_cosPsi_l"), v0.mLambda(), v0.pt(), angle * TMath::Cos(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosTheta_sinPhi_l"), v0.mLambda(), v0.pt(), angle * TMath::Sin(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_sinPsi_l"), v0.mLambda(), v0.pt(), angle * TMath::Sin(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosPhi_sinPsi_l"), v0.mLambda(), v0.pt(), TMath::Cos(v0.phi() * 2.0) * TMath::Sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/sinPhi_cosPsi_l"), v0.mLambda(), v0.pt(), TMath::Sin(v0.phi() * 2.0) * TMath::Cos(psi * 2.0), centrality);
          }
          if (cfgCalcCum1) {
            histos.fill(HIST("psi2/QA/cosTheta_l"), v0.mLambda(), v0.pt(), angle, centrality);
            histos.fill(HIST("psi2/QA/cosPsi_l"), v0.mLambda(), v0.pt(), TMath::Cos(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosPhi_l"), v0.mLambda(), v0.pt(), TMath::Cos(v0.phi() * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosPhi_cosPsi_l"), v0.mLambda(), v0.pt(), TMath::Cos(v0.phi() * 2.0) * TMath::Cos(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_cosPhi_l"), v0.mLambda(), v0.pt(), angle * TMath::Cos(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_cosPsi_l"), v0.mLambda(), v0.pt(), angle * TMath::Cos(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/sinPhi_sinPsi_l"), v0.mLambda(), v0.pt(), TMath::Sin(v0.phi() * 2.0) * TMath::Sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_sinPhi_l"), v0.mLambda(), v0.pt(), angle * TMath::Sin(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_sinPsi_l"), v0.mLambda(), v0.pt(), angle * TMath::Sin(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/sinPsi_l"), v0.mLambda(), v0.pt(), TMath::Sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/sinPhi_l"), v0.mLambda(), v0.pt(), TMath::Sin(v0.phi() * 2.0), centrality);
          }
        }
        if (aLambdaTag) {
          histos.fill(HIST("psi2/h_alambda_cos"), v0.mAntiLambda(), v0.pt(), angle * weight, centrality, relphi);
          histos.fill(HIST("psi2/h_alambda_cos2"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi2/h_alambda_cossin"), v0.mAntiLambda(), v0.pt(), angle * TMath::Sin(relphi) * weight, centrality);
          histos.fill(HIST("psi2/h_alambda_vncos"), v0.mAntiLambda(), v0.pt(), qvecMag * TMath::Cos(relphi) * weight, centrality);
          histos.fill(HIST("psi2/h_alambda_vnsin"), v0.mAntiLambda(), v0.pt(), TMath::Sin(relphi), centrality);

          if (cfgRapidityDep) {
            histos.fill(HIST("psi2/h_alambda_cos2_rap"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, v0.yLambda(), weight);
          }

          if (cfgAccAzimuth) {
            histos.fill(HIST("psi2/h_alambda_coscos"), v0.mAntiLambda(), v0.pt(), angle * TMath::Cos(relphi), centrality, weight);
          }

          if (cfgCalcCum) {
            histos.fill(HIST("psi2/QA/cosTheta_al"), v0.mAntiLambda(), v0.pt(), angle, centrality);
            histos.fill(HIST("psi2/QA/cosPsi_al"), v0.mAntiLambda(), v0.pt(), TMath::Cos(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosPhi_al"), v0.mAntiLambda(), v0.pt(), TMath::Cos(v0.phi() * 2.0), centrality);

            histos.fill(HIST("psi2/QA/sinPsi_al"), v0.mAntiLambda(), v0.pt(), TMath::Sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/sinPhi_al"), v0.mAntiLambda(), v0.pt(), TMath::Sin(v0.phi() * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosTheta_cosPhi_al"), v0.mAntiLambda(), v0.pt(), angle * TMath::Cos(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_cosPsi_al"), v0.mAntiLambda(), v0.pt(), angle * TMath::Cos(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosTheta_sinPhi_al"), v0.mAntiLambda(), v0.pt(), angle * TMath::Sin(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_sinPsi_al"), v0.mAntiLambda(), v0.pt(), angle * TMath::Sin(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosPhi_sinPsi_al"), v0.mAntiLambda(), v0.pt(), TMath::Cos(v0.phi() * 2.0) * TMath::Sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/sinPhi_cosPsi_al"), v0.mAntiLambda(), v0.pt(), TMath::Sin(v0.phi() * 2.0) * TMath::Cos(psi * 2.0), centrality);
          }
          if (cfgCalcCum1) {
            histos.fill(HIST("psi2/QA/cosTheta_al"), v0.mAntiLambda(), v0.pt(), angle, centrality);
            histos.fill(HIST("psi2/QA/cosPsi_al"), v0.mAntiLambda(), v0.pt(), TMath::Cos(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosPhi_al"), v0.mAntiLambda(), v0.pt(), TMath::Cos(v0.phi() * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosPhi_cosPsi_al"), v0.mAntiLambda(), v0.pt(), TMath::Cos(v0.phi() * 2.0) * TMath::Cos(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_cosPhi_al"), v0.mAntiLambda(), v0.pt(), angle * TMath::Cos(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_cosPsi_al"), v0.mAntiLambda(), v0.pt(), angle * TMath::Cos(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/sinPhi_sinPsi_al"), v0.mAntiLambda(), v0.pt(), TMath::Sin(v0.phi() * 2.0) * TMath::Sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_sinPhi_al"), v0.mAntiLambda(), v0.pt(), angle * TMath::Sin(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_sinPsi_al"), v0.mAntiLambda(), v0.pt(), angle * TMath::Sin(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/sinPsi_al"), v0.mAntiLambda(), v0.pt(), TMath::Sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/sinPhi_al"), v0.mAntiLambda(), v0.pt(), TMath::Sin(v0.phi() * 2.0), centrality);
          }
        }
      } else if (nmode == 3) {
        if (LambdaTag) {
          histos.fill(HIST("psi3/h_lambda_cos"), v0.mLambda(), v0.pt(), angle * weight, centrality, relphi);
          histos.fill(HIST("psi3/h_lambda_cos2"), v0.mLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi3/h_lambda_cossin"), v0.mLambda(), v0.pt(), angle * TMath::Sin(relphi) * weight, centrality);
          histos.fill(HIST("psi3/h_lambda_vncos"), v0.mLambda(), v0.pt(), qvecMag * TMath::Cos(relphi) * weight, centrality);
          histos.fill(HIST("psi3/h_lambda_vnsin"), v0.mLambda(), v0.pt(), TMath::Sin(relphi), centrality);

          if (cfgRapidityDep) {
            histos.fill(HIST("psi3/h_lambda_cos2_rap"), v0.mLambda(), v0.pt(), angle * angle, centrality, v0.yLambda(), weight);
          }

          if (cfgAccAzimuth) {
            histos.fill(HIST("psi3/h_lambda_coscos"), v0.mLambda(), v0.pt(), angle * TMath::Cos(relphi), centrality, weight);
          }
        }
        if (aLambdaTag) {
          histos.fill(HIST("psi3/h_alambda_cos"), v0.mAntiLambda(), v0.pt(), angle * weight, centrality, relphi);
          histos.fill(HIST("psi3/h_alambda_cos2"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, relphi, weight);
          histos.fill(HIST("psi3/h_alambda_cossin"), v0.mAntiLambda(), v0.pt(), angle * TMath::Sin(relphi) * weight, centrality);
          histos.fill(HIST("psi3/h_alambda_vncos"), v0.mAntiLambda(), v0.pt(), qvecMag * TMath::Cos(relphi) * weight, centrality);
          histos.fill(HIST("psi3/h_alambda_vnsin"), v0.mAntiLambda(), v0.pt(), TMath::Sin(relphi), centrality);

          if (cfgRapidityDep) {
            histos.fill(HIST("psi3/h_alambda_cos2_rap"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, v0.yLambda(), weight);
          }

          if (cfgAccAzimuth) {
            histos.fill(HIST("psi3/h_alambda_coscos"), v0.mAntiLambda(), v0.pt(), angle * TMath::Cos(relphi), centrality, weight);
          }
        }
      } else if (nmode == 4) {
        if (LambdaTag) {
          histos.fill(HIST("psi4/h_lambda_cos"), v0.mLambda(), v0.pt(), angle * weight, centrality, relphi);
          histos.fill(HIST("psi4/h_lambda_cos2"), v0.mLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi4/h_lambda_cossin"), v0.mLambda(), v0.pt(), angle * TMath::Sin(relphi) * weight, centrality);
          histos.fill(HIST("psi4/h_lambda_vncos"), v0.mLambda(), v0.pt(), qvecMag * TMath::Cos(relphi) * weight, centrality);
          histos.fill(HIST("psi4/h_lambda_vnsin"), v0.mLambda(), v0.pt(), TMath::Sin(relphi), centrality);

          if (cfgRapidityDep) {
            histos.fill(HIST("psi4/h_lambda_cos2_rap"), v0.mLambda(), v0.pt(), angle * angle, centrality, v0.yLambda(), weight);
          }

          if (cfgAccAzimuth) {
            histos.fill(HIST("psi4/h_lambda_coscos"), v0.mLambda(), v0.pt(), angle * TMath::Cos(relphi), centrality, weight);
          }
        }
        if (aLambdaTag) {
          histos.fill(HIST("psi4/h_alambda_cos"), v0.mAntiLambda(), v0.pt(), angle * weight, centrality, relphi);
          histos.fill(HIST("psi4/h_alambda_cos2"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi4/h_alambda_cossin"), v0.mAntiLambda(), v0.pt(), angle * TMath::Sin(relphi) * weight, centrality);
          histos.fill(HIST("psi4/h_alambda_vncos"), v0.mAntiLambda(), v0.pt(), qvecMag * TMath::Cos(relphi) * weight, centrality);
          histos.fill(HIST("psi4/h_alambda_vnsin"), v0.mAntiLambda(), v0.pt(), TMath::Sin(relphi), centrality);

          if (cfgRapidityDep) {
            histos.fill(HIST("psi4/h_alambda_cos2_rap"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, v0.yLambda(), weight);
          }

          if (cfgAccAzimuth) {
            histos.fill(HIST("psi4/h_alambda_coscos"), v0.mAntiLambda(), v0.pt(), angle * TMath::Cos(relphi), centrality, weight);
          }
        }
      } ////////// FIXME: not possible to get histograms using nmode
    }
  }

  void processData(EventCandidates::iterator const& collision,
                   TrackCandidates const& /*tracks*/, aod::V0Datas const& V0s,
                   aod::BCsWithTimestamps const&)
  {
    if (cfgCentEst == 1) {
      centrality = collision.centFT0C();
    } else if (cfgCentEst == 2) {
      centrality = collision.centFT0M();
    }
    if (!eventSelected(collision)) {
      return;
    }
    histos.fill(HIST("QA/CentDist"), centrality, 1.0);
    histos.fill(HIST("QA/PVzDist"), collision.posZ(), 1.0);

    if (cfgShiftCorr) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      currentRunNumber = bc.runNumber();
      if (currentRunNumber != lastRunNumber) {
        shiftprofile.clear();
        for (int i = 2; i < cfgnMods + 2; i++) {
          fullCCDBShiftCorrPath = cfgShiftPath;
          fullCCDBShiftCorrPath += "/v";
          fullCCDBShiftCorrPath += std::to_string(i);
          auto objshift = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPath, bc.timestamp());
          shiftprofile.push_back(objshift);
        }
        lastRunNumber = currentRunNumber;
      }
    }
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (cfgEffCor) {
      EffMap = ccdb->getForTimeStamp<TProfile2D>(cfgEffCorPath.value, bc.timestamp());
    }
    if (cfgAccCor) {
      AccMap = ccdb->getForTimeStamp<TProfile2D>(cfgAccCorPath.value, bc.timestamp());
    }
    for (int i = 2; i < cfgnMods + 2; i++) {
      if (cfgShiftCorrDef) {
        FillShiftCorrection(collision, i);
      }
      if (cfgQAv0) {
        FillEPQA(collision, i);
      }
      FillHistograms(collision, V0s, i);
    } // FIXME: need to fill different histograms for different harmonic
  }
  PROCESS_SWITCH(lambdapolarization, processData, "Process Event for data", true);

  using recoTracks = soa::Join<aod::TracksIU, aod::TracksExtra>;
  void processMC_ITSTPC(aod::McCollision const& mcCollision, soa::Join<aod::McParticles, aod::ParticlesToTracks> const& mcParticles, recoTracks const&)
  {
    float imp = mcCollision.impactParameter();
    float evPhi = mcCollision.eventPlaneAngle() / 2.0;
    float centclass = -999;
    if (imp >= 0 && imp < 3.49) {
      centclass = 2.5;
    }
    if (imp >= 3.49 && imp < 4.93) {
      centclass = 7.5;
    }
    if (imp >= 4.93 && imp < 6.98) {
      centclass = 15.0;
    }
    if (imp >= 6.98 && imp < 8.55) {
      centclass = 25.0;
    }
    if (imp >= 8.55 && imp < 9.87) {
      centclass = 35.0;
    }
    if (imp >= 9.87 && imp < 11) {
      centclass = 45.0;
    }
    if (imp >= 11 && imp < 12.1) {
      centclass = 55.0;
    }
    if (imp >= 12.1 && imp < 13.1) {
      centclass = 65.0;
    }
    if (imp >= 13.1 && imp < 14) {
      centclass = 75.0;
    }
    // if (evPhi < 0)
    //   evPhi += 2. * TMath::Pi();

    int nCh = 0;

    if (centclass > 0 && centclass < 80) {
      // event within range
      histos.fill(HIST("hImpactParameter"), imp);
      histos.fill(HIST("hEventPlaneAngle"), evPhi);
      for (auto const& mcParticle : mcParticles) {

        float deltaPhi = mcParticle.phi() - mcCollision.eventPlaneAngle();
        // focus on bulk: e, mu, pi, k, p
        int pdgCode = TMath::Abs(mcParticle.pdgCode());
        if (pdgCode != 3122)
          continue;
        if (!mcParticle.isPhysicalPrimary())
          continue;
        if (TMath::Abs(mcParticle.eta()) > 0.8) // main acceptance
          continue;
        histos.fill(HIST("hSparseMCGenWeight"), centclass, GetPhiInRange(deltaPhi), TMath::Power(TMath::Cos(2.0 * GetPhiInRange(deltaPhi)), 2.0), mcParticle.pt(), mcParticle.eta());
        nCh++;
        bool validGlobal = false;
        bool validAny = false;
        if (mcParticle.has_tracks()) {
          auto const& tracks = mcParticle.tracks_as<recoTracks>();
          for (auto const& track : tracks) {
            if (track.hasTPC() && track.hasITS()) {
              validGlobal = true;
            }
            if (track.hasTPC() || track.hasITS()) {
              validAny = true;
            }
          }
        }
        // if valid global, fill
        if (validGlobal) {
          histos.fill(HIST("hSparseMCRecWeight"), centclass, GetPhiInRange(deltaPhi), TMath::Power(TMath::Cos(2.0 * GetPhiInRange(deltaPhi)), 2.0), mcParticle.pt(), mcParticle.eta());
        }
        if (validAny) {
          histos.fill(HIST("hSparseMCRecAllTrackWeight"), centclass, GetPhiInRange(deltaPhi), TMath::Power(TMath::Cos(2.0 * GetPhiInRange(deltaPhi)), 2.0), mcParticle.pt(), mcParticle.eta());
          histos.fill(HIST("hEventPlaneAngleRec"), GetPhiInRange(deltaPhi));
        }
        // if any track present, fill
      }
    }
    histos.fill(HIST("hNchVsImpactParameter"), imp, nCh);
  }
  PROCESS_SWITCH(lambdapolarization, processMC_ITSTPC, "Process MC for ITSTPC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdapolarization>(cfgc, TaskName{"lf-lambdapolarization"})};
}
