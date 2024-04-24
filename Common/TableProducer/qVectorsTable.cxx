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
/// \file   qVectorsTable.cxx
/// \author Cindy Mordasini <cindy.mordasini@cern.ch>
/// \author Anna Önnerstad <anna.onnerstad@cern.ch>
///
/// \brief  Task calculating the Q-vectors for each collision in a bunch crossing
///         (with or without corrections) and save the results in a dedicated table.
///

// C++/ROOT includes.
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>

// o2Physics includes.
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"

#include "Common/DataModel/Qvectors.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
// o2 includes.
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

using namespace o2;
using namespace o2::framework;

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected,
                               aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As>;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension>;

struct qVectorsTable {
  enum {
    kFT0C = 0,
    kFT0A = 1,
    kFT0M,
    kFV0A,
    kBPos,
    kBNeg
  };

  // Configurables.
  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL",
                                     "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than",
                                      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(),
                                      "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  Configurable<int> cfgCentEsti{"cfgCentEsti",
                                2, "Centrality estimator (Run3): 0 = FT0M, 1 = FT0A, 2 = FT0C, 3 = FV0A"};

  Configurable<int> cfgCCDBConst{"cfgCCDBConst", 1, "Using constants in CCDB, 1 = CCDB, 2= Configurable"};
  Configurable<int> cfgGainCor{"cfgGainCor", 2, "Gain equalization, 1 = CCDB, 2 = Configurable"};
  // LOKI: We have here all centrality estimators for Run 3 (except FDDM and NTPV),
  // but the Q-vectors are calculated only for some of them.
  // FIXME: 6 correction factors for each centrality and 8 centrality intervals are hard-coded.
  // TODO: Constants from the CCDB
  Configurable<std::vector<float>> cfgFT0CCorr{"cfgFT0CCorr", std::vector<float>{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0C"};
  Configurable<std::vector<float>> cfgFT0ACorr{"cfgFT0ACorr", std::vector<float>{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0A"};
  Configurable<std::vector<float>> cfgFT0MCorr{"cfgFT0MCorr", std::vector<float>{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FT0M"};
  Configurable<std::vector<float>> cfgFV0ACorr{"cfgFV0ACorr", std::vector<float>{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for FV0A"};
  Configurable<std::vector<float>> cfgBPosCorr{"cfgBPosCorr", std::vector<float>{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for positive TPC tracks"};
  Configurable<std::vector<float>> cfgBNegCorr{"cfgBNegCorr", std::vector<float>{0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, "Correction constants for negative TPC tracks"};

  Configurable<std::vector<float>> cfgFT0RelGain{"cfgFT0RelGain", std::vector<float>{1.66815, 1.07726, 1.06665, 0.804114, 0.583846, 0.500518, 0.391628, 0.348394, 0.358615, 0.308253, 0.294405, 0.273764, 2.57031, 2.39421, 1.44468, 1.36256, 0.680442, 0.514505, 0.402177, 0.0768064, 0.724138, 0.69854, 0.45071, 0.468864, 2.39613, 1.29564, 1.32493, 0.966852, 0.847806, 0.676617, 0.560133, 0.445884, 0.449895, 0.361821, 0.361393, 0.352685, 3.34893, 3.14722, 1.77013, 1.6832, 0.880161, 0.671814, 0.574672, 0.500735, 0.911163, 0.869385, 0.568519, 0.575432, 2.6472, 1.48855, 1.54706, 1.14217, 0.787196, 0.598142, 0.53321, 0.44489, 0.490051, 0.385222, 0.41518, 0.366924, 2.9316, 2.8146, 1.52534, 1.61405, 0.899687, 0.701258, 0.54537, 0.506217, 0.823043, 0.904671, 0.548924, 0.54579, 2.42676, 1.45846, 1.48897, 1.02953, 0.827955, 0.640462, 0.572353, 0.46783, 0.488863, 0.369599, 0.415494, 0.362218, 3.17981, 3.01309, 1.79391, 1.65753, 0.922038, 0.747622, 0.585332, 0.516699, 1.04287, 1.00833, 0.673634, 0.647385, 1.28287, 0.982116, 0.952414, 0.812895, 0.696049, 0.643981, 0.561084, 0.545641, 0.627786, 0.556424, 0.580068, 0.563328, 1.55845, 1.35077, 1.08229, 0.932524, 0.721666, 0.673458, 0.544954, 0.57362, 0.633485, 0.627168, 0.545195, 0.614894, 1.71862, 1.4596, 1.13659, 1.0249, 0.941048, 0.69596, 0.621792, 0.609313, 0.727359, 0.618647, 0.651608, 0.668898, 1.8986, 1.74193, 1.33445, 1.08025, 0.823063, 0.773975, 0.665728, 0.661659, 0.71767, 0.682773, 0.678768, 0.703515, 2.09321, 1.70391, 1.31288, 1.13727, 0.842259, 0.782933, 0.691555, 0.66877, 0.729401, 0.657522, 0.677497, 0.652054, 1.6339, 1.73831, 1.58303, 1.17792, 0.888915, 0.833191, 0.693254, 0.689346, 0.80103, 0.751452, 0.741275, 0.757127, 2.32919, 1.93853, 1.46963, 1.27367, 0.957618, 1.07039, 0.737812, 0.759759, 0.827746, 0.724172, 0.782507, 0.803106, 2.80548, 0.99413, 1.73022, 1.50227, 0.921537, 0.869511, 1.03225, 1.07005, 1.57744, 1.30007, 1.23155, 1.06504, 1.70968, 1.25775, 1.24086, 1.07188, 1.7137, 1.36342, 1.30506, 1.12737, 1.82987, 1.39909, 1.14134, 1, 1, 1, 1, 1, 1}, "constants for relative FT0 gain equalization"};
  Configurable<std::vector<float>> cfgFV0RelGain{"cfgFV0RelGain", std::vector<float>{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "constants for relative FV0 gain equalization"};

  Configurable<float> cfgMinPtOnTPC{"cfgMinPtOnTPC", 0.15, "minimum transverse momentum selection for TPC tracks participating in Q-vector reconstruction"};
  Configurable<float> cfgMaxPtOnTPC{"cfgMaxPtOnTPC", 5., "maximum transverse momentum selection for TPC tracks participating in Q-vector reconstruction"};
  Configurable<int> cfgnMod{"cfgnMod", 2, "Modulation of interest"};

  Configurable<std::string> cfgGainEqPath{"cfgGainEqPath", "Users/j/junlee/Qvector/GainEq", "CCDB path for gain equalization constants"};
  Configurable<std::string> cfgQvecCalibPath{"cfgQvecCalibPath", "Analysis/EventPlane/QVecCorrections", "CCDB pasth for Q-vecteor calibration constants"};

  ConfigurableAxis cfgaxisFITamp{"cfgaxisFITamp", {1000, 0, 5000}, ""};

  // Table.
  Produces<aod::Qvectors> qVector;
  Produces<aod::QvectorFT0Cs> qVectorFT0C;
  Produces<aod::QvectorFT0As> qVectorFT0A;
  Produces<aod::QvectorFT0Ms> qVectorFT0M;
  Produces<aod::QvectorFV0As> qVectorFV0A;
  Produces<aod::QvectorBPoss> qVectorBPos;
  Produces<aod::QvectorBNegs> qVectorBNeg;

  std::vector<std::vector<float>> cfgCorr{};
  std::vector<float> FT0RelGainConst{};
  std::vector<float> FV0RelGainConst{};

  // Enable access to the CCDB for the offset and correction constants and save them
  // in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // geometry instances for V0 and T0
  o2::fv0::Geometry* fv0geom;
  o2::ft0::Geometry ft0geom;

  // Variables for other classes.
  EventPlaneHelper helperEP;

  HistogramRegistry histosQA{"histosQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int runNumber{-1};

  void init(InitContext const&)
  {
    // Setup the access to the CCDB objects of interest.
    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(cfgCcdbParam.nolaterthan.value);
    ccdb->setFatalWhenNull(false);

    AxisSpec axisPt = {40, 0.0, 4.0};
    AxisSpec axisEta = {32, -0.8, 0.8};
    AxisSpec axisPhi = {32, 0, constants::math::TwoPI};
    AxisSpec axixCent = {20, 0, 100};

    AxisSpec axisFITamp{cfgaxisFITamp, "FIT amp"};
    AxisSpec axisChID = {220, 0, 220};

    fv0geom = o2::fv0::Geometry::instance(o2::fv0::Geometry::eUninitialized);

    histosQA.add("ChTracks", "", {HistType::kTHnSparseF, {axisPt, axisEta, axisPhi, axixCent}});
    histosQA.add("FT0Amp", "", {HistType::kTH2F, {axisFITamp, axisChID}});
    histosQA.add("FT0AmpCor", "", {HistType::kTH2F, {axisFITamp, axisChID}});
    histosQA.add("FV0Amp", "", {HistType::kTH2F, {axisFITamp, axisChID}});
    histosQA.add("FV0AmpCor", "", {HistType::kTH2F, {axisFITamp, axisChID}});
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    cfgCorr.clear();
    FT0RelGainConst.clear();
    FV0RelGainConst.clear();
    cfgCorr = {};
    FT0RelGainConst = {};
    FV0RelGainConst = {};

    auto timestamp = bc.timestamp();

    auto offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", timestamp);
    auto offsetFV0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FV0/Calib/Align", timestamp);

    if (offsetFT0 != nullptr) {
      helperEP.SetOffsetFT0A((*offsetFT0)[0].getX(), (*offsetFT0)[0].getY());
      helperEP.SetOffsetFT0C((*offsetFT0)[1].getX(), (*offsetFT0)[1].getY());
    } else {
      LOGF(fatal, "Could not get the alignment parameters for FT0.");
    }

    if (offsetFV0 != nullptr) {
      helperEP.SetOffsetFV0left((*offsetFV0)[0].getX(), (*offsetFV0)[0].getY());
      helperEP.SetOffsetFV0right((*offsetFV0)[1].getX(), (*offsetFV0)[1].getY());
    } else {
      LOGF(fatal, "Could not get the alignment parameters for FV0.");
    }

    std::string fullPath;
    if (cfgCCDBConst == 1) {
      fullPath = cfgQvecCalibPath;
      fullPath += "/FT0C";
      auto objft0c = ccdb->getForTimeStamp<std::vector<float>>(fullPath, timestamp);
      if (!objft0c) {
        if (cfgFT0CCorr->size() < 48) {
          LOGF(fatal, "No proper correction factor assigned for FT0C");
        } else {
          cfgCorr.push_back(cfgFT0CCorr);
        }
      } else {
        cfgCorr.push_back(*(objft0c));
      }

      fullPath = cfgQvecCalibPath;
      fullPath += "/FT0A";
      auto objft0a = ccdb->getForTimeStamp<std::vector<float>>(fullPath, timestamp);
      if (!objft0a) {
        if (cfgFT0ACorr->size() < 48) {
          LOGF(fatal, "No proper correction factor assigned for FT0A");
        } else {
          cfgCorr.push_back(cfgFT0ACorr);
        }
      } else {
        cfgCorr.push_back(*(objft0a));
      }

      fullPath = cfgQvecCalibPath;
      fullPath += "/FT0M";
      auto objft0m = ccdb->getForTimeStamp<std::vector<float>>(fullPath, timestamp);
      if (!objft0m) {
        if (cfgFT0MCorr->size() < 48) {
          LOGF(fatal, "No proper correction factor assigned for FT0M");
        } else {
          cfgCorr.push_back(cfgFT0MCorr);
        }
      } else {
        cfgCorr.push_back(*(objft0m));
      }

      fullPath = cfgQvecCalibPath;
      fullPath += "/FV0A";
      auto objfv0a = ccdb->getForTimeStamp<std::vector<float>>(fullPath, timestamp);
      if (!objfv0a) {
        if (cfgFV0ACorr->size() < 48) {
          LOGF(fatal, "No proper correction factor assigned for FV0A");
        } else {
          cfgCorr.push_back(cfgFV0ACorr);
        }
      } else {
        cfgCorr.push_back(*(objfv0a));
      }

      fullPath = cfgQvecCalibPath;
      fullPath += "/BPos";
      auto objbpos = ccdb->getForTimeStamp<std::vector<float>>(fullPath, timestamp);
      if (!objbpos) {
        if (cfgBPosCorr->size() < 48) {
          LOGF(fatal, "No proper correction factor assigned for BPos");
        } else {
          cfgCorr.push_back(cfgBPosCorr);
        }
      } else {
        cfgCorr.push_back(*(objbpos));
      }

      fullPath = cfgQvecCalibPath;
      fullPath += "/BNeg";
      auto objbneg = ccdb->getForTimeStamp<std::vector<float>>(fullPath, timestamp);
      if (!objbneg) {
        if (cfgBNegCorr->size() < 48) {
          LOGF(fatal, "No proper correction factor assigned for BNeg");
        } else {
          cfgCorr.push_back(cfgBNegCorr);
        }
      } else {
        cfgCorr.push_back(*(objbneg));
      }
    } else if (cfgCCDBConst == 2) {
      if (cfgFT0CCorr->size() < 48) {
        LOGF(fatal, "No proper correction factor assigned for FT0C");
      }
      if (cfgFT0ACorr->size() < 48) {
        LOGF(fatal, "No proper correction factor assigned for FT0A");
      }
      if (cfgFT0MCorr->size() < 48) {
        LOGF(fatal, "No proper correction factor assigned for FT0M");
      }
      if (cfgFV0ACorr->size() < 48) {
        LOGF(fatal, "No proper correction factor assigned for FV0A");
      }
      if (cfgBPosCorr->size() < 48) {
        LOGF(fatal, "No proper correction factor assigned for positive TPC tracks");
      }
      if (cfgBNegCorr->size() < 48) {
        LOGF(fatal, "No proper correction factor assigned for negative TPC tracks");
      } // will be replaced with method that call constants from CCDB

      cfgCorr.push_back(cfgFT0CCorr);
      cfgCorr.push_back(cfgFT0ACorr);
      cfgCorr.push_back(cfgFT0MCorr);
      cfgCorr.push_back(cfgFV0ACorr);
      cfgCorr.push_back(cfgBPosCorr);
      cfgCorr.push_back(cfgBNegCorr);
    }

    if (cfgGainCor == 0) {
      for (auto i{0u}; i < cfgFT0RelGain->size(); i++) {
        FT0RelGainConst.push_back(1.);
      }
      for (auto i{0u}; i < cfgFV0RelGain->size(); i++) {
        FV0RelGainConst.push_back(1.);
      }
    } else if (cfgGainCor == 1) {
      fullPath = cfgGainEqPath;
      fullPath += "/FT0";
      auto objft0Gain = ccdb->getForTimeStamp<std::vector<float>>(fullPath, timestamp);
      if (!objft0Gain) {
        for (auto i{0u}; i < cfgFT0RelGain->size(); i++) {
          FT0RelGainConst.push_back(1.);
        }
      } else {
        FT0RelGainConst = *(objft0Gain);
      }

      fullPath = cfgGainEqPath;
      fullPath += "/FV0";
      auto objfv0Gain = ccdb->getForTimeStamp<std::vector<float>>(fullPath, timestamp);
      if (!objfv0Gain) {
        for (auto i{0u}; i < cfgFV0RelGain->size(); i++) {
          FV0RelGainConst.push_back(1.);
        }
      } else {
        FV0RelGainConst = *(objfv0Gain);
      }
    } else if (cfgGainCor == 2) {
      FT0RelGainConst = cfgFT0RelGain;
      FV0RelGainConst = cfgFV0RelGain;
    }
  }

  template <typename TrackType>
  bool SelTrack(const TrackType track)
  {
    if (track.pt() < cfgMinPtOnTPC)
      return false;
    if (track.pt() > cfgMaxPtOnTPC)
      return false;
    if (!track.passedITSNCls())
      return false;
    if (!track.passedITSChi2NDF())
      return false;
    if (!track.passedITSHits())
      return false;
    if (!track.passedTPCCrossedRowsOverNCls())
      return false;
    if (!track.passedTPCChi2NDF())
      return false;
    if (!track.passedDCAxy())
      return false;
    if (!track.passedDCAz())
      return false;

    return true;
  }

  void process(MyCollisions::iterator const& coll, aod::BCsWithTimestamps const&, aod::FT0s const&, aod::FV0As const&, MyTracks const& tracks)
  {

    std::vector<int> TrkBPosLabel{};
    std::vector<int> TrkBNegLabel{};
    std::vector<float> qvecRe{};
    std::vector<float> qvecIm{};
    std::vector<float> qvecAmp{};

    auto bc = coll.bc_as<aod::BCsWithTimestamps>();
    int currentRun = bc.runNumber();
    if (runNumber != currentRun) {
      initCCDB(bc);
      runNumber = currentRun;
    }

    // Get the centrality value for all subscribed estimators and takes the one
    // corresponding to cfgCentEsti. Reject also the events with invalid values.
    // NOTE: centFDDM and centNTPV not implemented as it makes the compilation crashes...
    float centAllEstim[4] = {
      coll.centFT0M(), coll.centFT0A(), coll.centFT0C(),
      coll.centFV0A()};
    float cent = centAllEstim[cfgCentEsti];
    if (cent < 0. || cent > 100.) {
      cent = 110.;
    }

    // Calculate the Q-vectors values for this event.
    // TODO: Add here qVect for other detectors,...
    float qVectFT0A[2] = {0.}; // Real and imaginary parts of the Q-vector in FT0A.
    float qVectFT0C[2] = {0.}; // Real and imaginary parts of the Q-vector in FT0C.
    float qVectFT0M[2] = {0.};
    float qVectFV0A[2] = {0.}; // Real and imaginary parts of the Q-vector in FV0A.

    float qVectBPos[2] = {0.};
    float qVectBNeg[2] = {0.};

    TComplex QvecDet(0); // Complex value of the Q-vector for any detector.
    TComplex QvecFT0M(0);
    float sumAmplFT0A = 0.; // Sum of the amplitudes of all non-dead channels in any detector.
    float sumAmplFT0C = 0.;
    float sumAmplFT0M = 0.;
    float sumAmplFV0A = 0.;

    /// First check if the collision has a found FT0. If yes, calculate the
    /// Q-vectors for FT0A and FT0C (both real and imaginary parts). If no,
    /// attribute dummy values to the corresponding qVect.
    if (coll.has_foundFT0()) {
      auto ft0 = coll.foundFT0();

      // Iterate over the non-dead channels for FT0-A to get the total Q-vector
      // and sum of amplitudes.
      for (std::size_t iChA = 0; iChA < ft0.channelA().size(); iChA++) {
        // Get first the corresponding amplitude.
        float ampl = ft0.amplitudeA()[iChA];
        int FT0AchId = ft0.channelA()[iChA];

        histosQA.fill(HIST("FT0Amp"), ampl, FT0AchId);
        histosQA.fill(HIST("FT0AmpCor"), ampl / FT0RelGainConst[FT0AchId], FT0AchId);
        // Update the Q-vector and sum of amplitudes using the helper function.
        // LOKI: Note this assumes nHarmo = 2!! Likely generalise in the future.
        helperEP.SumQvectors(0, FT0AchId, ampl / FT0RelGainConst[FT0AchId], cfgnMod, QvecDet, sumAmplFT0A, ft0geom, fv0geom);
        helperEP.SumQvectors(0, FT0AchId, ampl / FT0RelGainConst[FT0AchId], cfgnMod, QvecFT0M, sumAmplFT0M, ft0geom, fv0geom);
      } // Go to the next channel iChA.

      // Set the Qvectors for FT0A with the normalised Q-vector values if the sum of
      // amplitudes is non-zero. Otherwise, set it to a dummy 999.
      if (sumAmplFT0A > 1e-8) {
        QvecDet /= sumAmplFT0A;
        qVectFT0A[0] = QvecDet.Re();
        qVectFT0A[1] = QvecDet.Im();
        // printf("qVectFT0A[0] = %.2f ; qVectFT0A[1] = %.2f \n", qVectFT0A[0], qVectFT0A[1]); // Debug printing.
      } else {
        qVectFT0A[0] = 999.;
        qVectFT0A[1] = 999.;
      }

      // Repeat the procedure with FT0-C for the found FT0.
      // Start by resetting to zero the intermediate quantities.
      QvecDet = TComplex(0., 0.);
      for (std::size_t iChC = 0; iChC < ft0.channelC().size(); iChC++) {
        // iChC ranging from 0 to max 112. We need to add 96 (= max channels in FT0-A)
        // to ensure a proper channel number in FT0 as a whole.
        float ampl = ft0.amplitudeC()[iChC];
        int FT0CchId = ft0.channelC()[iChC] + 96;

        histosQA.fill(HIST("FT0Amp"), ampl, FT0CchId);
        histosQA.fill(HIST("FT0AmpCor"), ampl / FT0RelGainConst[FT0CchId], FT0CchId);

        helperEP.SumQvectors(0, FT0CchId, ampl / FT0RelGainConst[FT0CchId], cfgnMod, QvecDet, sumAmplFT0C, ft0geom, fv0geom);
        helperEP.SumQvectors(0, FT0CchId, ampl / FT0RelGainConst[FT0CchId], cfgnMod, QvecFT0M, sumAmplFT0M, ft0geom, fv0geom);
      }

      if (sumAmplFT0C > 1e-8) {
        QvecDet /= sumAmplFT0C;
        qVectFT0C[0] = QvecDet.Re();
        qVectFT0C[1] = QvecDet.Im();
        // printf("qVectFT0C[0] = %.2f ; qVectFT0C[1] = %.2f \n", qVectFT0C[0], qVectFT0C[1]); // Debug printing.
      } else {
        qVectFT0C[0] = 999.;
        qVectFT0C[1] = 999.;
      }

      if (sumAmplFT0M > 1e-8) {
        QvecFT0M /= sumAmplFT0M;
        qVectFT0M[0] = QvecFT0M.Re();
        qVectFT0M[1] = QvecFT0M.Im();
      } else {
        qVectFT0M[0] = 999.;
        qVectFT0M[1] = 999.;
      }
    } else {
      qVectFT0A[0] = -999.;
      qVectFT0A[1] = -999.;
      qVectFT0C[0] = -999.;
      qVectFT0C[1] = -999.;
      qVectFT0M[0] = -999.;
      qVectFT0M[1] = -999.;
    }

    QvecDet = TComplex(0., 0.);
    sumAmplFV0A = 0;
    if (coll.has_foundFV0()) {
      auto fv0 = coll.foundFV0();

      for (std::size_t iCh = 0; iCh < fv0.channel().size(); iCh++) {
        float ampl = fv0.amplitude()[iCh];
        int FV0AchId = fv0.channel()[iCh];
        histosQA.fill(HIST("FV0Amp"), ampl, FV0AchId);
        histosQA.fill(HIST("FV0AmpCor"), ampl / FV0RelGainConst[FV0AchId], FV0AchId);

        helperEP.SumQvectors(1, FV0AchId, ampl / FV0RelGainConst[FV0AchId], cfgnMod, QvecDet, sumAmplFV0A, ft0geom, fv0geom);
      }

      if (sumAmplFV0A > 1e-8) {
        QvecDet /= sumAmplFV0A;
        qVectFV0A[0] = QvecDet.Re();
        qVectFV0A[1] = QvecDet.Im();
        // printf("qVectFV0[0] = %.2f ; qVectFV0[1] = %.2f \n", qVectFV0[0], qVectFV0[1]); // Debug printing.
      } else {
        qVectFV0A[0] = 999.;
        qVectFV0A[1] = 999.;
      }
    } else {
      qVectFV0A[0] = -999.;
      qVectFV0A[1] = -999.;
    }

    int nTrkBPos = 0;
    int nTrkBNeg = 0;

    for (auto& trk : tracks) {
      if (!SelTrack(trk)) {
        continue;
      }
      histosQA.fill(HIST("ChTracks"), trk.pt(), trk.eta(), trk.phi(), cent);
      if (std::abs(trk.eta()) < 0.1 || std::abs(trk.eta()) > 0.8) {
        continue;
      }
      if (trk.eta() > 0) {
        qVectBPos[0] += trk.pt() * std::cos(trk.phi() * cfgnMod);
        qVectBPos[1] += trk.pt() * std::sin(trk.phi() * cfgnMod);
        TrkBPosLabel.push_back(trk.globalIndex());
        nTrkBPos++;
      } else if (trk.eta() < 0) {
        qVectBNeg[0] += trk.pt() * std::cos(trk.phi() * cfgnMod);
        qVectBNeg[1] += trk.pt() * std::sin(trk.phi() * cfgnMod);
        TrkBNegLabel.push_back(trk.globalIndex());
        nTrkBNeg++;
      }
    }

    if (nTrkBPos > 0) {
      qVectBPos[0] /= nTrkBPos;
      qVectBPos[1] /= nTrkBPos;
    } else {
      qVectBPos[0] = 999.;
      qVectBPos[1] = 999.;
    }

    if (nTrkBNeg > 0) {
      qVectBNeg[0] /= nTrkBNeg;
      qVectBNeg[1] /= nTrkBNeg;
    } else {
      qVectBNeg[0] = 999.;
      qVectBNeg[1] = 999.;
    }

    int cBin = helperEP.GetCentBin(cent);

    for (auto i{0u}; i < 4; i++) {
      qvecRe.push_back(qVectFT0C[0]);
      qvecIm.push_back(qVectFT0C[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      qvecRe.push_back(qVectFT0A[0]);
      qvecIm.push_back(qVectFT0A[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      qvecRe.push_back(qVectFT0M[0]);
      qvecIm.push_back(qVectFT0M[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      qvecRe.push_back(qVectFV0A[0]);
      qvecIm.push_back(qVectFV0A[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      qvecRe.push_back(qVectBPos[0]);
      qvecIm.push_back(qVectBPos[1]);
    }
    for (auto i{0u}; i < 4; i++) {
      qvecRe.push_back(qVectBNeg[0]);
      qvecIm.push_back(qVectBNeg[1]);
    }

    qvecAmp.push_back(sumAmplFT0C);
    qvecAmp.push_back(sumAmplFT0A);
    qvecAmp.push_back(sumAmplFT0M);
    qvecAmp.push_back(sumAmplFV0A);
    qvecAmp.push_back(static_cast<float>(nTrkBPos));
    qvecAmp.push_back(static_cast<float>(nTrkBNeg));

    if (cBin != -1) {
      for (auto i{0u}; i < 6; i++) {
        helperEP.DoRecenter(qvecRe[i * 4 + 1], qvecIm[i * 4 + 1], cfgCorr[i][cBin * 6], cfgCorr[i][cBin * 6 + 1]);

        helperEP.DoRecenter(qvecRe[i * 4 + 2], qvecIm[i * 4 + 2], cfgCorr[i][cBin * 6], cfgCorr[i][cBin * 6 + 1]);
        helperEP.DoTwist(qvecRe[i * 4 + 2], qvecIm[i * 4 + 2], cfgCorr[i][cBin * 6 + 2], cfgCorr[i][cBin * 6 + 3]);

        helperEP.DoRecenter(qvecRe[i * 4 + 3], qvecIm[i * 4 + 3], cfgCorr[i][cBin * 6], cfgCorr[i][cBin * 6 + 1]);
        helperEP.DoTwist(qvecRe[i * 4 + 3], qvecIm[i * 4 + 3], cfgCorr[i][cBin * 6 + 2], cfgCorr[i][cBin * 6 + 3]);
        helperEP.DoRescale(qvecRe[i * 4 + 3], qvecIm[i * 4 + 3], cfgCorr[i][cBin * 6 + 4], cfgCorr[i][cBin * 6 + 5]);
      }
    }

    // Fill the columns of the Qvectors table.
    qVector(cent, cBin, qvecRe, qvecIm, qvecAmp);
    qVectorFT0C(cBin, qvecRe[kFT0C * 4 + 3], qvecIm[kFT0C * 4 + 3], sumAmplFT0C);
    qVectorFT0A(cBin, qvecRe[kFT0A * 4 + 3], qvecIm[kFT0A * 4 + 3], sumAmplFT0A);
    qVectorFT0M(cBin, qvecRe[kFT0M * 4 + 3], qvecIm[kFT0M * 4 + 3], sumAmplFT0M);
    qVectorFV0A(cBin, qvecRe[kFV0A * 4 + 3], qvecIm[kFV0A * 4 + 3], sumAmplFV0A);
    qVectorBPos(cBin, qvecRe[kBPos * 4 + 3], qvecIm[kBPos * 4 + 3], nTrkBPos, TrkBPosLabel);
    qVectorBNeg(cBin, qvecRe[kBNeg * 4 + 3], qvecIm[kBNeg * 4 + 3], nTrkBNeg, TrkBNegLabel);

  } // End process.
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qVectorsTable>(cfgc)};
}
