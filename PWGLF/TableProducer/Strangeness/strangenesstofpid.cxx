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
//
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//    Strangeness TOF PID
//  *+-+*+-+*+-+*+-+*+-+*+-+*
//
/// \author Nicolò Jacazio
/// \author David Dobrigkeit Chinellato
/// \since  11/05/2023
/// \brief  Table producer for V0 daughter PID info
//
// This task produces daughter PID information for strange daughters
// taking into account the (candidate-by-candidate) time spent as a heavier
// (strange, weakly-decaying) particle. This task is meant to be a test, as
// it hasn't been fully tested yet! Use at your own peril for now :-)

#include "TableHelper.h"

#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <map>
#include <utility>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// For original data loops
using V0OriginalDatas = soa::Join<aod::V0Indices, aod::V0Cores>;
using CascOriginalDatas = soa::Join<aod::CascIndices, aod::CascCores>;
using TracksWithAllExtras = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::TOFEvTime, aod::TOFSignal>;

// For derived data analysis
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using V0DerivedDatas = soa::Join<aod::V0Cores, aod::V0Extras, aod::V0CollRefs>;
using V0DerivedDatasMC = soa::Join<aod::V0Cores, aod::V0Extras, aod::V0CollRefs, aod::V0CoreMCLabels>;
using CascDerivedDatas = soa::Join<aod::CascCores, aod::CascExtras, aod::CascCollRefs>;
using CascDerivedDatasMC = soa::Join<aod::CascCores, aod::CascExtras, aod::CascCollRefs, aod::CascCoreMCLabels>;

struct strangenesstofpid {
  // TOF pid for strangeness (recalculated with topology)
  Produces<aod::V0TOFPIDs> v0tofpid;            // table with Nsigmas
  Produces<aod::V0TOFBetas> v0tofbeta;          // table with betas
  Produces<aod::V0TOFDebugs> v0tofdebugs;       // table with extra debug information
  Produces<aod::V0TOFNSigmas> v0tofnsigmas;     // table with nsigmas
  Produces<aod::CascTOFPIDs> casctofpids;       // cascades: table with base info
  Produces<aod::CascTOFNSigmas> casctofnsigmas; // cascades: table with Nsigmas

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // mean vertex position to be used if no collision associated
  o2::dataformats::MeanVertexObject* mVtx = nullptr;

  // LUT for Propagator + TrackLTIntegral
  o2::base::MatLayerCylSet* lut = nullptr;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // master switches
  Configurable<int> calculationMethod{"calculationMethod", 0, "algorithm for TOF calculation. 0: fast analytical withouot eloss, 1: O2 Propagator + trackLTIntegral (slow), 2: both methods and do comparison studies (slow)"};
  Configurable<int> calculateV0s{"calculateV0s", -1, "calculate V0-related TOF PID (0: no, 1: yes, -1: auto)"};
  Configurable<int> calculateCascades{"calculateCascades", -1, "calculate cascade-related TOF PID (0: no, 1: yes, -1: auto)"};

  // Operation and minimisation criteria
  struct : ConfigurableGroup {
    Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
    Configurable<float> tofPosition{"tofPosition", 377.934f, "TOF effective (inscribed) radius"};
  } propagationConfiguration;

  Configurable<bool> doQA{"doQA", false, "create QA histos"};
  Configurable<bool> doNSigmas{"doNSigmas", true, "calculate TOF N-sigma"};
  Configurable<bool> doQANSigma{"doQANSigma", false, "create QA of Nsigma histos"};

  // configurables related to V0s
  struct : ConfigurableGroup {
    std::string prefix = "v0Calibration";
    Configurable<float> qaDCADau{"qaDCADau", 0.5, "DCA daughters (cm) for QA plots"};
    Configurable<float> qaCosPA{"qaCosPA", 0.999, "CosPA for QA plots"};
    Configurable<float> qaMassWindow{"qaMassWindow", 0.005, "Mass window around expected (in GeV/c2) for QA plots"};
    Configurable<float> qaTPCNSigma{"qaTPCNSigma", 5, "TPC N-sigma to apply for qa plots"};
  } v0Group;

  // configurables related to V0s
  struct : ConfigurableGroup {
    std::string prefix = "cascadeCalibration";
    Configurable<float> qaV0DCADau{"qaV0DCADau", 0.5, "DCA daughters (cm) for QA plots"};
    Configurable<float> qaCascDCADau{"qaCascDCADau", 0.5, "DCA daughters (cm) for QA plots"};
    Configurable<float> qaV0CosPA{"qaV0CosPA", 0.995, "CosPA for QA plots"};
    Configurable<float> qaCascCosPA{"qaCascCosPA", 0.995, "CosPA for QA plots"};
    Configurable<float> qaMassWindow{"qaMassWindow", 0.005, "Mass window around expected (in GeV/c2) for QA plots"};
    Configurable<float> qaTPCNSigma{"qaTPCNSigma", 5, "TPC N-sigma to apply for qa plots"};
  } cascadeGroup;

  // CCDB options
  // CCDB options
  struct : ConfigurableGroup {
    std::string prefix = "ccdb";
    Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> nSigmaPath{"nSigmaPath", "Users/d/ddobrigk/stratof", "Path of information for n-sigma calculation"};
    Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};
  } ccdbConfigurations;

  // manual
  Configurable<int> useCustomRunNumber{"useCustomRunNumber", false, "Use custom timestamp"};
  Configurable<int> manualRunNumber{"manualRunNumber", 544122, "manual run number if no collisions saved"};

  ConfigurableAxis axisPosition{"axisPosition", {400, -400.f, +400.f}, "position (cm)"};
  ConfigurableAxis axisEta{"axisEta", {20, -1.0f, +1.0f}, "#eta"};
  ConfigurableAxis axisDeltaTime{"axisDeltaTime", {2000, -1000.0f, +1000.0f}, "delta-time (ps)"};
  ConfigurableAxis axisTime{"axisTime", {400, 10000.0f, +50000.0f}, "T (ps)"};
  ConfigurableAxis axisNSigma{"axisNSigma", {200, -10.0f, +10.0f}, "N(#sigma)"};
  ConfigurableAxis axisRatioMethods{"axisRatioMethods", {400, 0.9f, 1.9f}, "T_{method 1}/T_{method 0}"};
  ConfigurableAxis axisSnp{"axisSnp", {220, -1.1f, 1.1f}, "snp"};

  // master p axis
  ConfigurableAxis axisP{"axisP", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};

  // for zooming in at low values only (e-loss studies and effective correction)
  ConfigurableAxis axisSmallP{"axisSmallP", {250, 0.0f, 2.5f}, "p_{T} (GeV/c)"};

  // for n-sigma calibration
  bool nSigmaCalibLoaded;
  TList* nSigmaCalibObjects = nullptr;
  TH1 *hMeanPosLaPi = nullptr, *hSigmaPosLaPi = nullptr;
  TH1 *hMeanPosLaPr = nullptr, *hSigmaPosLaPr = nullptr;
  TH1 *hMeanNegLaPi = nullptr, *hSigmaNegLaPi = nullptr;
  TH1 *hMeanNegLaPr = nullptr, *hSigmaNegLaPr = nullptr;
  TH1 *hMeanPosK0Pi = nullptr, *hSigmaPosK0Pi = nullptr;
  TH1 *hMeanNegK0Pi = nullptr, *hSigmaNegK0Pi = nullptr;
  TH1 *hMeanPosXiPi = nullptr, *hSigmaPosXiPi = nullptr;
  TH1 *hMeanPosXiPr = nullptr, *hSigmaPosXiPr = nullptr;
  TH1 *hMeanNegXiPi = nullptr, *hSigmaNegXiPi = nullptr;
  TH1 *hMeanNegXiPr = nullptr, *hSigmaNegXiPr = nullptr;
  TH1 *hMeanBachXiPi = nullptr, *hSigmaBachXiPi = nullptr;
  TH1 *hMeanPosOmPi = nullptr, *hSigmaPosOmPi = nullptr;
  TH1 *hMeanPosOmPr = nullptr, *hSigmaPosOmPr = nullptr;
  TH1 *hMeanNegOmPi = nullptr, *hSigmaNegOmPi = nullptr;
  TH1 *hMeanNegOmPr = nullptr, *hSigmaNegOmPr = nullptr;
  TH1 *hMeanBachOmKa = nullptr, *hSigmaBachOmKa = nullptr;

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation

  // enum to keep track of the TOF-related properties for V0s
  enum tofEnum { kLength = 0,
                 kHasTOF,
                 kNEnums };

  /// function to calculate track length of this track up to a certain segment of a detector
  /// to be used internally in another function that calculates length until it finds the proper one
  /// warning: this could be optimised further for speed
  /// \param track the input track
  /// \param x1 x of the first point of the detector segment
  /// \param y1 y of the first point of the detector segment
  /// \param x2 x of the first point of the detector segment
  /// \param y2 y of the first point of the detector segment
  /// \param magneticField the magnetic field to use when propagating
  float trackLengthToSegment(o2::track::TrackPar track, float x1, float y1, float x2, float y2, float magneticField)
  {
    // don't make use of the track parametrization
    float length = -104;

    // causality protection
    std::array<float, 3> mom;
    track.getPxPyPzGlo(mom);
    // get start point
    std::array<float, 3> startPoint;
    track.getXYZGlo(startPoint);

    // better replaced with scalar momentum check later
    // if (((x1 + x2) * mom[0] + (y1 + y2) * mom[1]) < 0.0f)
    //   return -101;

    // get circle X, Y please
    o2::math_utils::CircleXYf_t trcCircle;
    float sna, csa;
    track.getCircleParams(magneticField, trcCircle, sna, csa);

    // Calculate necessary inner product
    float segmentModulus = std::hypot(x2 - x1, y2 - y1);
    float alongSegment = ((trcCircle.xC - x1) * (x2 - x1) + (trcCircle.yC - y1) * (y2 - y1)) / segmentModulus;

    // find point of closest approach between segment and circle center
    float pcaX = (x2 - x1) * alongSegment / segmentModulus + x1;
    float pcaY = (y2 - y1) * alongSegment / segmentModulus + y1;

    float centerDistToPC = std::hypot(pcaX - trcCircle.xC, pcaY - trcCircle.yC);

    // distance pca-to-intercept in multiples of segment modulus (for convenience)
    if (centerDistToPC > trcCircle.rC)
      return -103;

    float pcaToIntercept = TMath::Sqrt(TMath::Abs(trcCircle.rC * trcCircle.rC - centerDistToPC * centerDistToPC));

    float interceptX1 = pcaX + (x2 - x1) / segmentModulus * pcaToIntercept;
    float interceptY1 = pcaY + (y2 - y1) / segmentModulus * pcaToIntercept;
    float interceptX2 = pcaX - (x2 - x1) / segmentModulus * pcaToIntercept;
    float interceptY2 = pcaY - (y2 - y1) / segmentModulus * pcaToIntercept;

    float scalarCheck1 = ((x2 - x1) * (interceptX1 - x1) + (y2 - y1) * (interceptY1 - y1)) / segmentModulus;
    float scalarCheck2 = ((x2 - x1) * (interceptX2 - x1) + (y2 - y1) * (interceptY2 - y1)) / segmentModulus;

    float cosAngle1 = -1000, sinAngle1 = -1000, modulus1 = -1000;
    float cosAngle2 = -1000, sinAngle2 = -1000, modulus2 = -1000;
    float length1 = -1000, length2 = -1000;

    modulus1 = std::hypot(interceptX1 - trcCircle.xC, interceptY1 - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
    cosAngle1 = (interceptX1 - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (interceptY1 - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
    sinAngle1 = (interceptX1 - trcCircle.xC) * (startPoint[1] - trcCircle.yC) - (interceptY1 - trcCircle.yC) * (startPoint[0] - trcCircle.xC);
    cosAngle1 /= modulus1;
    sinAngle1 /= modulus1;
    length1 = trcCircle.rC * TMath::ACos(cosAngle1);
    length1 *= sqrt(1.0f + track.getTgl() * track.getTgl());

    modulus2 = std::hypot(interceptX2 - trcCircle.xC, interceptY2 - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
    cosAngle2 = (interceptX2 - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (interceptY2 - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
    sinAngle2 = (interceptX2 - trcCircle.xC) * (startPoint[1] - trcCircle.yC) - (interceptY2 - trcCircle.yC) * (startPoint[0] - trcCircle.xC);
    cosAngle2 /= modulus2;
    sinAngle2 /= modulus2;
    length2 = trcCircle.rC * TMath::ACos(cosAngle2);
    length2 *= sqrt(1.0f + track.getTgl() * track.getTgl());

    // rotate transverse momentum vector such that it is at intercepts
    float angle1 = TMath::ACos(cosAngle1);
    if (sinAngle1 < 0)
      angle1 *= -1.0f;
    float px1 = +TMath::Cos(angle1) * mom[0] + TMath::Sin(angle1) * mom[1];
    float py1 = -TMath::Sin(angle1) * mom[0] + TMath::Cos(angle1) * mom[1];

    float angle2 = TMath::ACos(cosAngle2);
    if (sinAngle2 < 0)
      angle2 *= -1.0f;
    float px2 = +TMath::Cos(angle2) * mom[0] + TMath::Sin(angle2) * mom[1];
    float py2 = -TMath::Sin(angle2) * mom[0] + TMath::Cos(angle2) * mom[1];

    float midSegX = 0.5f * (x2 + x1);
    float midSegY = 0.5f * (y2 + y1);

    float scalarMomentumCheck1 = px1 * midSegX + py1 * midSegY;
    float scalarMomentumCheck2 = px2 * midSegX + py2 * midSegY;

    if (scalarCheck1 > 0.0f && scalarCheck1 < segmentModulus && scalarMomentumCheck1 > 0.0f) {
      length = length1;
      // X = interceptX1; Y = interceptY1;
    }
    if (scalarCheck2 > 0.0f && scalarCheck2 < segmentModulus && scalarMomentumCheck2 > 0.0f) {
      length = length2;
      // X = interceptX2; Y = interceptY2;
    }
    return length;
  }

  /// function to calculate track length of this track up to a certain segmented detector
  /// \param track the input track
  /// \param magneticField the magnetic field to use when propagating
  float findInterceptLength(o2::track::TrackPar track, float magneticField)
  {
    float length = 1e+6;
    for (int iSeg = 0; iSeg < 18; iSeg++) {
      // Detector segmentation loop
      float segmentAngle = 20.0f / 180.0f * TMath::Pi();
      float theta = static_cast<float>(iSeg) * 20.0f / 180.0f * TMath::Pi();
      float halfWidth = propagationConfiguration.tofPosition * TMath::Tan(0.5f * segmentAngle);
      float x1 = TMath::Cos(theta) * (-halfWidth) + TMath::Sin(theta) * propagationConfiguration.tofPosition;
      float y1 = -TMath::Sin(theta) * (-halfWidth) + TMath::Cos(theta) * propagationConfiguration.tofPosition;
      float x2 = TMath::Cos(theta) * (+halfWidth) + TMath::Sin(theta) * propagationConfiguration.tofPosition;
      float y2 = -TMath::Sin(theta) * (+halfWidth) + TMath::Cos(theta) * propagationConfiguration.tofPosition;
      float thisLength = trackLengthToSegment(track, x1, y1, x2, y2, magneticField);
      if (thisLength < length && thisLength > 0) {
        length = thisLength;
      }
    }
    if (length > 1e+5)
      length = -100; // force negative to avoid misunderstandings
    return length;
  }

  void init(InitContext& initContext)
  {
    if (calculateV0s.value < 0) {
      // check if TOF information is required, enable if so
      calculateV0s.value = isTableRequiredInWorkflow(initContext, "V0TOFNSigmas");
      if (calculateV0s.value > 0) {
        LOGF(info, "Strangeness TOF PID: V0 calculations enabled automatically");
      }
    }
    if (calculateCascades.value < 0) {
      // check if TOF information is required, enable if so
      calculateCascades.value = isTableRequiredInWorkflow(initContext, "CascTOFNSigmas");
      if (calculateCascades.value > 0) {
        LOGF(info, "Strangeness TOF PID: Cascade calculations enabled automatically");
      }
    }

    nSigmaCalibLoaded = false;
    nSigmaCalibObjects = nullptr;

    // for n-sigma calibration
    hMeanPosLaPi = nullptr;
    hSigmaPosLaPi = nullptr;
    hMeanPosLaPr = nullptr;
    hSigmaPosLaPr = nullptr;
    hMeanNegLaPi = nullptr;
    hSigmaNegLaPi = nullptr;
    hMeanNegLaPr = nullptr;
    hSigmaNegLaPr = nullptr;
    hMeanPosK0Pi = nullptr;
    hSigmaNegK0Pi = nullptr;
    hMeanNegK0Pi = nullptr;
    hSigmaNegK0Pi = nullptr;

    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  // could be changed later
    maxStep = 2.00f; // could be changed later

    ccdb->setURL(ccdbConfigurations.ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // per event
    histos.add("hCandidateCounter", "hCandidateCounter", kTH1F, {{500, -0.5f, 499.5f}});

    // measured vs expected total time QA
    if (doQA) {
      // standard deltaTime values
      if (calculateV0s.value > 0) {
        histos.add("h2dDeltaTimePositiveLambdaPi", "h2dDeltaTimePositiveLambdaPi", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dDeltaTimeNegativeLambdaPi", "h2dDeltaTimeNegativeLambdaPi", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dDeltaTimePositiveLambdaPr", "h2dDeltaTimePositiveLambdaPr", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dDeltaTimeNegativeLambdaPr", "h2dDeltaTimeNegativeLambdaPr", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dDeltaTimePositiveK0ShortPi", "h2dDeltaTimePositiveK0ShortPi", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dDeltaTimeNegativeK0ShortPi", "h2dDeltaTimeNegativeK0ShortPi", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
      }

      if (calculateCascades.value > 0) {
        histos.add("h2dposDeltaTimeAsXiPi", "h2dposDeltaTimeAsXiPi", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dposDeltaTimeAsXiPr", "h2dposDeltaTimeAsXiPr", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dnegDeltaTimeAsXiPi", "h2dnegDeltaTimeAsXiPi", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dnegDeltaTimeAsXiPr", "h2dnegDeltaTimeAsXiPr", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dbachDeltaTimeAsXiPi", "h2dbachDeltaTimeAsXiPi", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});

        histos.add("h2dposDeltaTimeAsOmPi", "h2dposDeltaTimeAsOmPi", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dposDeltaTimeAsOmPr", "h2dposDeltaTimeAsOmPr", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dnegDeltaTimeAsOmPi", "h2dnegDeltaTimeAsOmPi", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dnegDeltaTimeAsOmPr", "h2dnegDeltaTimeAsOmPr", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
        histos.add("h2dbachDeltaTimeAsOmKa", "h2dbachDeltaTimeAsOmKa", {HistType::kTH3F, {axisP, axisEta, axisDeltaTime}});
      }

      histos.add("h2dPositiveTOFProperties", "h2dPositiveTOFProperties", {HistType::kTH2F, {axisP, {4, -0.5, 3.5f}}});
      histos.add("h2dNegativeTOFProperties", "h2dNegativeTOFProperties", {HistType::kTH2F, {axisP, {4, -0.5, 3.5f}}});

      if (doQANSigma) {
        if (calculateV0s.value > 0) {
          histos.add("h2dNSigmaPositiveLambdaPi", "h2dNSigmaPositiveLambdaPi", {HistType::kTH2F, {axisP, axisNSigma}});
          histos.add("h2dNSigmaNegativeLambdaPi", "h2dNSigmaNegativeLambdaPi", {HistType::kTH2F, {axisP, axisNSigma}});
          histos.add("h2dNSigmaPositiveLambdaPr", "h2dNSigmaPositiveLambdaPr", {HistType::kTH2F, {axisP, axisNSigma}});
          histos.add("h2dNSigmaNegativeLambdaPr", "h2dNSigmaNegativeLambdaPr", {HistType::kTH2F, {axisP, axisNSigma}});
          histos.add("h2dNSigmaPositiveK0ShortPi", "h2dNSigmaPositiveK0ShortPi", {HistType::kTH2F, {axisP, axisNSigma}});
          histos.add("h2dNSigmaNegativeK0ShortPi", "h2dNSigmaNegativeK0ShortPi", {HistType::kTH2F, {axisP, axisNSigma}});
        }

        if (calculateCascades.value > 0) {
          histos.add("h2dNSigmaXiLaPi", "h2dNSigmaXiLaPi", {HistType::kTH2F, {axisP, axisNSigma}});
          histos.add("h2dNSigmaXiLaPr", "h2dNSigmaXiLaPr", {HistType::kTH2F, {axisP, axisNSigma}});
          histos.add("h2dNSigmaXiPi", "h2dNSigmaXiPi", {HistType::kTH2F, {axisP, axisNSigma}});
          histos.add("h2dNSigmaOmLaPi", "h2dNSigmaOmLaPi", {HistType::kTH2F, {axisP, axisNSigma}});
          histos.add("h2dNSigmaOmLaPr", "h2dNSigmaOmLaPr", {HistType::kTH2F, {axisP, axisNSigma}});
          histos.add("h2dNSigmaOmKa", "h2dNSigmaOmKa", {HistType::kTH2F, {axisP, axisNSigma}});
        }
      }

      // delta lambda decay time
      histos.add("h2dLambdaDeltaDecayTime", "h2dLambdaDeltaDecayTime", {HistType::kTH2F, {axisP, axisDeltaTime}});
    }
  }

  void initCCDB(int runNumber)
  {
    if (mRunNumber == runNumber) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (propagationConfiguration.d_bz_input > -990) {
      d_bz = propagationConfiguration.d_bz_input;
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mVtx = ccdb->getForRun<o2::dataformats::MeanVertexObject>(ccdbConfigurations.mVtxPath, runNumber);
      mRunNumber = runNumber;
      return;
    }

    o2::parameters::GRPObject* grpo = ccdb->getForRun<o2::parameters::GRPObject>(ccdbConfigurations.grpPath, runNumber);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for run " << runNumber << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(ccdbConfigurations.grpmagPath, runNumber);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << ccdbConfigurations.grpmagPath << " of object GRPMagField and " << ccdbConfigurations.grpPath << " of object GRPObject for run " << runNumber;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      mVtx = ccdb->getForRun<o2::dataformats::MeanVertexObject>(ccdbConfigurations.mVtxPath, runNumber);
      LOG(info) << "Retrieved GRP for run " << runNumber << " with magnetic field of " << d_bz << " kZG";
    }

    // if TOF Nsigma desired
    if (doNSigmas) {
      nSigmaCalibObjects = ccdb->getForRun<TList>(ccdbConfigurations.nSigmaPath, runNumber);
      if (nSigmaCalibObjects) {
        LOGF(info, "loaded TList with this many objects: %i", nSigmaCalibObjects->GetEntries());
        nSigmaCalibLoaded = true; // made it thus far, mark loaded

        if (calculateV0s.value) {
          hMeanPosLaPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanPosLaPi"));
          hMeanPosLaPr = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanPosLaPr"));
          hMeanNegLaPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanNegLaPi"));
          hMeanNegLaPr = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanNegLaPr"));
          hMeanPosK0Pi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanPosK0Pi"));
          hMeanNegK0Pi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanNegK0Pi"));

          hSigmaPosLaPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaPosLaPi"));
          hSigmaPosLaPr = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaPosLaPr"));
          hSigmaNegLaPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaNegLaPi"));
          hSigmaNegLaPr = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaNegLaPr"));
          hSigmaPosK0Pi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaPosK0Pi"));
          hSigmaNegK0Pi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaNegK0Pi"));

          if (!hMeanPosLaPi)
            LOG(info) << "Problems finding mean histogram hMeanPosLaPi!";
          if (!hMeanPosLaPr)
            LOG(info) << "Problems finding mean histogram hMeanPosLaPr!";
          if (!hMeanNegLaPi)
            LOG(info) << "Problems finding mean histogram hMeanNegLaPi!";
          if (!hMeanNegLaPr)
            LOG(info) << "Problems finding mean histogram hMeanNegLaPr!";
          if (!hMeanPosK0Pi)
            LOG(info) << "Problems finding mean histogram hMeanPosK0Pi!";
          if (!hMeanNegK0Pi)
            LOG(info) << "Problems finding mean histogram hMeanNegK0Pi!";
          if (!hSigmaPosK0Pi || !hSigmaNegK0Pi || !hSigmaPosLaPi || !hSigmaPosLaPr || !hSigmaNegLaPi || !hSigmaNegLaPr) {
            LOG(info) << "Problems finding sigma histograms!";
          }

          if (calculateCascades.value) {
            hMeanPosXiPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanPosXiPi"));
            hMeanPosXiPr = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanPosXiPr"));
            hMeanNegXiPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanNegXiPi"));
            hMeanNegXiPr = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanNegXiPr"));
            hMeanBachXiPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanBachXiPi"));
            hMeanPosOmPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanPosOmPi"));
            hMeanPosOmPr = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanPosOmPr"));
            hMeanNegOmPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanNegOmPi"));
            hMeanNegOmPr = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanNegOmPr"));
            hMeanBachOmKa = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hMeanBachOmKa"));

            hSigmaPosXiPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaPosXiPi"));
            hSigmaPosXiPr = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaPosXiPr"));
            hSigmaNegXiPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaNegXiPi"));
            hSigmaNegXiPr = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaNegXiPr"));
            hSigmaBachXiPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaBachXiPi"));
            hSigmaPosOmPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaPosOmPi"));
            hSigmaPosOmPr = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaPosOmPr"));
            hSigmaNegOmPi = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaNegOmPi"));
            hSigmaNegOmPr = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaNegOmPr"));
            hSigmaBachOmKa = reinterpret_cast<TH1*>(nSigmaCalibObjects->FindObject("hSigmaBachOmKa"));

            if (!hMeanPosXiPi || !hMeanPosXiPr || !hMeanNegXiPi || !hMeanNegXiPr || !hMeanBachXiPi)
              LOG(info) << "Problems finding xi mean histograms!";
            if (!hMeanPosOmPi || !hMeanPosOmPr || !hMeanNegOmPi || !hMeanNegOmPr || !hMeanBachOmKa)
              LOG(info) << "Problems finding omega sigma histograms!";
            if (!hSigmaPosXiPi || !hSigmaPosXiPr || !hSigmaNegXiPi || !hSigmaNegXiPr || !hSigmaBachXiPi)
              LOG(info) << "Problems finding xi sigma histograms!";
            if (!hSigmaPosOmPi || !hSigmaPosOmPr || !hSigmaNegOmPi || !hSigmaNegOmPr || !hSigmaBachOmKa)
              LOG(info) << "Problems finding omega sigma histograms!";
          }
        }
      }
    }

    // if (calculationMethod.value > 0 && !lut) {
    //   // setMatLUT only after magfield has been initalized
    //   // (setMatLUT has implicit and problematic init field call if not)
    //   LOG(info) << "Loading full (all-radius) material look-up table for run number: " << runNumber;
    //   lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForRun<o2::base::MatLayerCylSet>(ccdbConfigurations.lutPath, runNumber));
    //   o2::base::Propagator::Instance()->setMatLUT(lut);
    //   o2::base::Propagator::Instance()->setTGeoFallBackAllowed(false);
    //   LOG(info) << "Material look-up table loaded!";
    // }
    mRunNumber = runNumber;
  }

  float velocity(float lMomentum, float lMass)
  {
    // Momentum p and mass m -> returns speed in centimeters per picosecond
    // Useful for TOF calculations
    float lA = (lMomentum / lMass) * (lMomentum / lMass);
    return 0.0299792458 * TMath::Sqrt(lA / (1 + lA));
  }

  // structs to hold information
  struct v0TofInfo { // holds processed information regarding TOF for V0s
    float timePositivePr = o2::aod::v0data::kNoTOFValue;
    float timePositivePi = o2::aod::v0data::kNoTOFValue;
    float timeNegativePr = o2::aod::v0data::kNoTOFValue;
    float timeNegativePi = o2::aod::v0data::kNoTOFValue;

    float deltaTimePositiveLambdaPi = o2::aod::v0data::kNoTOFValue;
    float deltaTimeNegativeLambdaPi = o2::aod::v0data::kNoTOFValue;
    float deltaTimePositiveLambdaPr = o2::aod::v0data::kNoTOFValue;
    float deltaTimeNegativeLambdaPr = o2::aod::v0data::kNoTOFValue;
    float deltaTimePositiveK0ShortPi = o2::aod::v0data::kNoTOFValue;
    float deltaTimeNegativeK0ShortPi = o2::aod::v0data::kNoTOFValue;

    float nSigmaPositiveLambdaPi = o2::aod::v0data::kNoTOFValue;
    float nSigmaPositiveLambdaPr = o2::aod::v0data::kNoTOFValue;
    float nSigmaNegativeLambdaPi = o2::aod::v0data::kNoTOFValue;
    float nSigmaNegativeLambdaPr = o2::aod::v0data::kNoTOFValue;
    float nSigmaPositiveK0ShortPi = o2::aod::v0data::kNoTOFValue;
    float nSigmaNegativeK0ShortPi = o2::aod::v0data::kNoTOFValue;
  };

  // structs to hold information
  struct cascTofInfo { // holds processed information regarding TOF for Cascades
    float posFlightPi = o2::aod::cascdata::kNoTOFValue;
    float posFlightPr = o2::aod::cascdata::kNoTOFValue;
    float negFlightPi = o2::aod::cascdata::kNoTOFValue;
    float negFlightPr = o2::aod::cascdata::kNoTOFValue;
    float bachFlightPi = o2::aod::cascdata::kNoTOFValue;
    float bachFlightKa = o2::aod::cascdata::kNoTOFValue;

    float posDeltaTimeAsXiPi = o2::aod::cascdata::kNoTOFValue, posDeltaTimeAsXiPr = o2::aod::cascdata::kNoTOFValue;
    float negDeltaTimeAsXiPi = o2::aod::cascdata::kNoTOFValue, negDeltaTimeAsXiPr = o2::aod::cascdata::kNoTOFValue;
    float bachDeltaTimeAsXiPi = o2::aod::cascdata::kNoTOFValue;
    float posDeltaTimeAsOmPi = o2::aod::cascdata::kNoTOFValue, posDeltaTimeAsOmPr = o2::aod::cascdata::kNoTOFValue;
    float negDeltaTimeAsOmPi = o2::aod::cascdata::kNoTOFValue, negDeltaTimeAsOmPr = o2::aod::cascdata::kNoTOFValue;
    float bachDeltaTimeAsOmKa = o2::aod::cascdata::kNoTOFValue;

    float nSigmaXiLaPr = o2::aod::cascdata::kNoTOFValue;
    float nSigmaXiLaPi = o2::aod::cascdata::kNoTOFValue;
    float nSigmaXiPi = o2::aod::cascdata::kNoTOFValue;
    float nSigmaOmLaPr = o2::aod::cascdata::kNoTOFValue;
    float nSigmaOmLaPi = o2::aod::cascdata::kNoTOFValue;
    float nSigmaOmKa = o2::aod::cascdata::kNoTOFValue;
  };

  struct trackTofInfo { // holds input track info
    bool hasITS = false;
    bool hasTPC = false;
    bool hasTOF = false;
    int collisionId = -1;
    float tofExpMom = 0.0f;
    float tofSignal = 0.0f;
    float tofEvTime = 0.0f;
    float length = 0.0f;

    // save TPC PID here for completeness too
    float tpcNSigmaPi = 0.0f;
    float tpcNSigmaKa = 0.0f;
    float tpcNSigmaPr = 0.0f;
  };

  // templatized process function for symmetric operation in derived and original AO2D
  /// \param collisions the collisions table (needed for de-referencing V0 and progns)
  /// \param v0 the V0 being processed
  /// \param pTof the TOF information for the positive track
  /// \param nTof the TOF information for the negative track
  template <class TCollisions, typename TV0, typename TTOFInfo>
  v0TofInfo calculateTofInfoV0(TCollisions const& collisions, int const& collisionId, TV0 const& v0, TTOFInfo const& pTof, TTOFInfo const& nTof)
  {
    v0TofInfo v0tof; // return this struct
    auto collision = collisions.rawIteratorAt(collisionId);

    //_____________________________________________________________________________________________
    // daughter tracks: initialize from V0 position and momenta
    o2::track::TrackPar posTrack = o2::track::TrackPar({v0.x(), v0.y(), v0.z()}, {v0.pxpos(), v0.pypos(), v0.pzpos()}, +1, false);
    o2::track::TrackPar negTrack = o2::track::TrackPar({v0.x(), v0.y(), v0.z()}, {v0.pxneg(), v0.pyneg(), v0.pzneg()}, -1, false);

    //_____________________________________________________________________________________________
    // time of V0 segment
    float lengthV0 = std::hypot(v0.x() - collision.posX(), v0.y() - collision.posY(), v0.z() - collision.posZ());
    float velocityK0Short = velocity(v0.p(), o2::constants::physics::MassKaonNeutral);
    float velocityLambda = velocity(v0.p(), o2::constants::physics::MassLambda);
    float timeK0Short = lengthV0 / velocityK0Short; // in picoseconds
    float timeLambda = lengthV0 / velocityLambda;   // in picoseconds

    //_____________________________________________________________________________________________
    // define simple checks
    bool passesQAcuts = (v0.v0cosPA() > v0Group.qaCosPA && v0.dcaV0daughters() < v0Group.qaDCADau);
    bool lambdaCandidate = std::abs(v0.mLambda() - o2::constants::physics::MassLambda) < v0Group.qaMassWindow &&
                           std::abs(pTof.tpcNSigmaPr) < v0Group.qaTPCNSigma &&
                           std::abs(nTof.tpcNSigmaPi) < v0Group.qaTPCNSigma;
    bool antiLambdaCandidate = std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda) < v0Group.qaMassWindow &&
                               std::abs(pTof.tpcNSigmaPi) < v0Group.qaTPCNSigma &&
                               std::abs(nTof.tpcNSigmaPr) < v0Group.qaTPCNSigma;
    bool k0ShortCandidate = std::abs(v0.mK0Short() - o2::constants::physics::MassKaonNeutral) < v0Group.qaMassWindow &&
                            std::abs(pTof.tpcNSigmaPi) < v0Group.qaTPCNSigma &&
                            std::abs(nTof.tpcNSigmaPi) < v0Group.qaTPCNSigma;

    //_____________________________________________________________________________________________
    // Actual calculation
    if (pTof.hasTOF && pTof.hasITS) {
      float velocityPositivePr, velocityPositivePi, lengthPositive;
      velocityPositivePr = velocityPositivePi = lengthPositive = o2::aod::v0data::kNoTOFValue;
      // method 0: legacy standalone without use of primary particle TOF
      if (calculationMethod.value == 0) {
        velocityPositivePr = velocity(posTrack.getP(), o2::constants::physics::MassProton);
        velocityPositivePi = velocity(posTrack.getP(), o2::constants::physics::MassPionCharged);
        lengthPositive = findInterceptLength(posTrack, d_bz);
        v0tof.timePositivePr = lengthPositive / velocityPositivePr;
        v0tof.timePositivePi = lengthPositive / velocityPositivePi;
      }
      // method 1: correct primary particle TOF information
      // length -> revise by removing travel length to primary vertex
      // expected momentum -> kept as is for now, could correct at second stage
      // use main method from TOF to calculate expected time
      if (calculationMethod.value == 1) {
        if (pTof.collisionId >= 0) {
          auto trackCollision = collisions.rawIteratorAt(pTof.collisionId);
          const o2::math_utils::Point3D<float> trackVertex{trackCollision.posX(), trackCollision.posY(), trackCollision.posZ()};
          o2::track::TrackLTIntegral ltIntegral;
          bool successPropag = o2::base::Propagator::Instance()->propagateToDCA(trackVertex, posTrack, d_bz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrNONE, nullptr, &ltIntegral);
          if (successPropag) {
            lengthPositive = pTof.length - ltIntegral.getL();
            v0tof.timePositivePr = o2::framework::pid::tof::MassToExpTime(pTof.tofExpMom, lengthPositive, o2::constants::physics::MassProton * o2::constants::physics::MassProton);
            v0tof.timePositivePi = o2::framework::pid::tof::MassToExpTime(pTof.tofExpMom, lengthPositive, o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);
          }
        }
      }
      if (lengthPositive > 0.0f) {
        v0tof.deltaTimePositiveLambdaPr = (pTof.tofSignal - pTof.tofEvTime) - (timeLambda + v0tof.timePositivePr);
        v0tof.deltaTimePositiveLambdaPi = (pTof.tofSignal - pTof.tofEvTime) - (timeLambda + v0tof.timePositivePi);
        v0tof.deltaTimePositiveK0ShortPi = (pTof.tofSignal - pTof.tofEvTime) - (timeK0Short + v0tof.timePositivePi);

        // de facto nsigma
        if (nSigmaCalibLoaded) {
          v0tof.nSigmaPositiveLambdaPi = (v0tof.deltaTimePositiveLambdaPi - hMeanPosLaPi->Interpolate(v0.p())) / hSigmaPosLaPi->Interpolate(v0.p());
          v0tof.nSigmaPositiveLambdaPr = (v0tof.deltaTimePositiveLambdaPr - hMeanPosLaPr->Interpolate(v0.p())) / hSigmaPosLaPr->Interpolate(v0.p());
          v0tof.nSigmaPositiveK0ShortPi = (v0tof.deltaTimePositiveK0ShortPi - hMeanPosK0Pi->Interpolate(v0.p())) / hSigmaPosK0Pi->Interpolate(v0.p());
        }

        // do QA histograms (calibration / QC)
        if (doQA) {
          if (passesQAcuts) {
            if (lambdaCandidate) {
              histos.fill(HIST("h2dDeltaTimePositiveLambdaPr"), v0.p(), v0.eta(), v0tof.deltaTimePositiveLambdaPr);
            }
            if (antiLambdaCandidate) {
              histos.fill(HIST("h2dDeltaTimePositiveLambdaPi"), v0.p(), v0.eta(), v0tof.deltaTimePositiveLambdaPi);
            }
            if (k0ShortCandidate) {
              histos.fill(HIST("h2dDeltaTimePositiveK0ShortPi"), v0.p(), v0.eta(), v0tof.deltaTimePositiveK0ShortPi);
            }
          }
        }
      }
    }
    if (nTof.hasTOF && nTof.hasITS) {
      float velocityNegativePr, velocityNegativePi, lengthNegative;
      velocityNegativePr = velocityNegativePi = lengthNegative = o2::aod::v0data::kNoTOFValue;
      // method 0: legacy standalone without use of primary particle TOF
      if (calculationMethod.value == 0) {
        velocityNegativePr = velocity(negTrack.getP(), o2::constants::physics::MassProton);
        velocityNegativePi = velocity(negTrack.getP(), o2::constants::physics::MassPionCharged);
        lengthNegative = findInterceptLength(negTrack, d_bz);
        v0tof.timeNegativePr = lengthNegative / velocityNegativePr;
        v0tof.timeNegativePi = lengthNegative / velocityNegativePi;
      }
      // method 1: correct primary particle TOF information
      // length -> revise by removing travel length to primary vertex
      // expected momentum -> kept as is for now, could correct at second stage
      // use main method from TOF to calculate expected time
      if (calculationMethod.value == 1) {
        if (nTof.collisionId >= 0) {
          auto trackCollision = collisions.rawIteratorAt(nTof.collisionId);
          const o2::math_utils::Point3D<float> trackVertex{trackCollision.posX(), trackCollision.posY(), trackCollision.posZ()};
          o2::track::TrackLTIntegral ltIntegral;
          bool successPropag = o2::base::Propagator::Instance()->propagateToDCA(trackVertex, negTrack, d_bz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrNONE, nullptr, &ltIntegral);
          if (successPropag) {
            lengthNegative = nTof.length - ltIntegral.getL();
            v0tof.timeNegativePr = o2::framework::pid::tof::MassToExpTime(nTof.tofExpMom, nTof.length - ltIntegral.getL(), o2::constants::physics::MassProton * o2::constants::physics::MassProton);
            v0tof.timeNegativePi = o2::framework::pid::tof::MassToExpTime(nTof.tofExpMom, nTof.length - ltIntegral.getL(), o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);
          }
        }
      }
      if (lengthNegative > 0.0f) {
        v0tof.deltaTimeNegativeLambdaPr = (nTof.tofSignal - nTof.tofEvTime) - (timeLambda + v0tof.timeNegativePr);
        v0tof.deltaTimeNegativeLambdaPi = (nTof.tofSignal - nTof.tofEvTime) - (timeLambda + v0tof.timeNegativePi);
        v0tof.deltaTimeNegativeK0ShortPi = (nTof.tofSignal - nTof.tofEvTime) - (timeK0Short + v0tof.timeNegativePi);

        // de facto nsigma
        if (nSigmaCalibLoaded) {
          v0tof.nSigmaNegativeLambdaPi = (v0tof.deltaTimeNegativeLambdaPi - hMeanNegLaPi->Interpolate(v0.p())) / hSigmaNegLaPi->Interpolate(v0.p());
          v0tof.nSigmaNegativeLambdaPr = (v0tof.deltaTimeNegativeLambdaPr - hMeanNegLaPr->Interpolate(v0.p())) / hSigmaNegLaPr->Interpolate(v0.p());
          v0tof.nSigmaNegativeK0ShortPi = (v0tof.deltaTimeNegativeK0ShortPi - hMeanNegK0Pi->Interpolate(v0.p())) / hSigmaNegK0Pi->Interpolate(v0.p());
        }

        // do QA histograms (calibration / QC)
        if (doQA) {
          if (passesQAcuts) {
            if (lambdaCandidate) {
              histos.fill(HIST("h2dDeltaTimeNegativeLambdaPi"), v0.p(), v0.eta(), v0tof.deltaTimeNegativeLambdaPi);
            }
            if (antiLambdaCandidate) {
              histos.fill(HIST("h2dDeltaTimeNegativeLambdaPr"), v0.p(), v0.eta(), v0tof.deltaTimeNegativeLambdaPr);
            }
            if (k0ShortCandidate) {
              histos.fill(HIST("h2dDeltaTimeNegativeK0ShortPi"), v0.p(), v0.eta(), v0tof.deltaTimeNegativeK0ShortPi);
            }
          }
        }
      }
    }

    return v0tof;
  } // end calculation altogether

  template <class TCollisions, typename TCascade, typename TTOFInfo>
  cascTofInfo calculateTofInfoCascade(TCollisions const& collisions, int const& collisionId, TCascade const& cascade, TTOFInfo const& pTof, TTOFInfo const& nTof, TTOFInfo const& bTof)
  {
    cascTofInfo casctof; // return this struct
    auto collision = collisions.rawIteratorAt(collisionId);

    //_____________________________________________________________________________________________
    // daughter tracks: initialize from V0 position and momenta
    o2::track::TrackPar posTrack = o2::track::TrackPar({cascade.xlambda(), cascade.ylambda(), cascade.zlambda()}, {cascade.pxpos(), cascade.pypos(), cascade.pzpos()}, +1, false);
    o2::track::TrackPar negTrack = o2::track::TrackPar({cascade.xlambda(), cascade.ylambda(), cascade.zlambda()}, {cascade.pxneg(), cascade.pyneg(), cascade.pzneg()}, -1, false);
    o2::track::TrackPar bachTrack = o2::track::TrackPar({cascade.x(), cascade.y(), cascade.z()}, {cascade.pxbach(), cascade.pybach(), cascade.pzbach()}, cascade.sign(), false);
    o2::track::TrackPar cascTrack = o2::track::TrackPar({cascade.x(), cascade.y(), cascade.z()}, {cascade.px(), cascade.py(), cascade.pz()}, cascade.sign(), false);

    //_____________________________________________________________________________________________
    // time of V0 segment
    float velocityXi = velocity(cascTrack.getP(), o2::constants::physics::MassXiMinus);
    float velocityOm = velocity(cascTrack.getP(), o2::constants::physics::MassOmegaMinus);
    float velocityLa = velocity(std::hypot(cascade.pxlambda(), cascade.pylambda(), cascade.pzlambda()), o2::constants::physics::MassLambda);
    float lengthV0 = std::hypot(cascade.xlambda() - cascade.x(), cascade.ylambda() - cascade.y(), cascade.zlambda() - cascade.z());
    float lengthCascade = o2::aod::cascdata::kNoTOFValue;

    // cascade length (N.B. could be simpler via trackLTIntegral, kept with legacy calculation)
    const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};
    bool successPropag = o2::base::Propagator::Instance()->propagateToDCA(collVtx, cascTrack, d_bz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrNONE);
    float d = -1.0f;
    float linearToPV = std::hypot(cascade.x() - collision.posX(), cascade.y() - collision.posY(), cascade.z() - collision.posZ());
    if (successPropag) {
      std::array<float, 3> cascCloseToPVPosition;
      cascTrack.getXYZGlo(cascCloseToPVPosition);
      o2::math_utils::CircleXYf_t trcCircleCascade;
      float sna, csa;
      cascTrack.getCircleParams(d_bz, trcCircleCascade, sna, csa);

      // calculate 2D distance between two points
      d = std::hypot(cascade.x() - cascCloseToPVPosition[0], cascade.y() - cascCloseToPVPosition[1]);
      // d3d = std::hypot(cascade.x() - cascCloseToPVPosition[0], cascade.y() - cascCloseToPVPosition[1], cascade.z() - cascCloseToPVPosition[2]); // cross-check variable
      float sinThetaOverTwo = d / (2.0f * trcCircleCascade.rC);
      lengthCascade = 2.0f * trcCircleCascade.rC * TMath::ASin(sinThetaOverTwo);
      lengthCascade *= sqrt(1.0f + cascTrack.getTgl() * cascTrack.getTgl());
    }

    if (!successPropag) {
      lengthCascade = linearToPV; // if propagation failed, use linear estimate (optional: actually do not define?)
    }

    // flight times of decaying particles
    float lambdaFlight = lengthV0 / velocityLa;
    float xiFlight = lengthCascade / velocityXi;
    float omFlight = lengthCascade / velocityOm;

    //_____________________________________________________________________________________________
    // define simple checks
    bool passesQAcuts = (cascade.dcaV0daughters() < cascadeGroup.qaV0DCADau && cascade.dcacascdaughters() < cascadeGroup.qaCascDCADau && cascade.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > cascadeGroup.qaV0CosPA && cascade.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > cascadeGroup.qaCascCosPA);
    bool xiMinusCandidate = cascade.sign() < 0 &&
                            std::abs(cascade.mXi() - o2::constants::physics::MassXiMinus) < cascadeGroup.qaMassWindow &&
                            std::abs(pTof.tpcNSigmaPr) < cascadeGroup.qaTPCNSigma &&
                            std::abs(nTof.tpcNSigmaPi) < cascadeGroup.qaTPCNSigma &&
                            std::abs(bTof.tpcNSigmaPi) < cascadeGroup.qaTPCNSigma;
    bool xiPlusCandidate = cascade.sign() > 0 &&
                           std::abs(cascade.mXi() - o2::constants::physics::MassXiMinus) < cascadeGroup.qaMassWindow &&
                           std::abs(pTof.tpcNSigmaPi) < cascadeGroup.qaTPCNSigma &&
                           std::abs(nTof.tpcNSigmaPr) < cascadeGroup.qaTPCNSigma &&
                           std::abs(bTof.tpcNSigmaPi) < cascadeGroup.qaTPCNSigma;
    bool omegaMinusCandidate = cascade.sign() < 0 &&
                               std::abs(cascade.mOmega() - o2::constants::physics::MassOmegaMinus) < cascadeGroup.qaMassWindow &&
                               std::abs(pTof.tpcNSigmaPr) < cascadeGroup.qaTPCNSigma &&
                               std::abs(nTof.tpcNSigmaPi) < cascadeGroup.qaTPCNSigma &&
                               std::abs(bTof.tpcNSigmaKa) < cascadeGroup.qaTPCNSigma;
    bool omegaPlusCandidate = cascade.sign() > 0 &&
                              std::abs(cascade.mOmega() - o2::constants::physics::MassOmegaMinus) < cascadeGroup.qaMassWindow &&
                              std::abs(pTof.tpcNSigmaPi) < cascadeGroup.qaTPCNSigma &&
                              std::abs(nTof.tpcNSigmaPr) < cascadeGroup.qaTPCNSigma &&
                              std::abs(bTof.tpcNSigmaKa) < cascadeGroup.qaTPCNSigma;

    //_____________________________________________________________________________________________
    // Actual calculation
    if (pTof.hasTOF && pTof.hasITS) {
      float velocityPositivePr, velocityPositivePi, lengthPositive;
      velocityPositivePr = velocityPositivePi = lengthPositive = o2::aod::v0data::kNoTOFValue;
      if (calculationMethod.value == 0) {
        velocityPositivePr = velocity(posTrack.getP(), o2::constants::physics::MassProton);
        velocityPositivePi = velocity(posTrack.getP(), o2::constants::physics::MassPionCharged);
        lengthPositive = findInterceptLength(posTrack, d_bz);
        casctof.posFlightPr = lengthPositive / velocityPositivePr;
        casctof.posFlightPi = lengthPositive / velocityPositivePi;
      }
      // method 1: correct primary particle TOF information
      // length -> revise by removing travel length to primary vertex
      // expected momentum -> kept as is for now, could correct at second stage
      // use main method from TOF to calculate expected time
      if (calculationMethod.value == 1) {
        if (pTof.collisionId >= 0) {
          auto trackCollision = collisions.rawIteratorAt(pTof.collisionId);
          const o2::math_utils::Point3D<float> trackVertex{trackCollision.posX(), trackCollision.posY(), trackCollision.posZ()};
          o2::track::TrackLTIntegral ltIntegral;
          bool successPropag = o2::base::Propagator::Instance()->propagateToDCA(trackVertex, posTrack, d_bz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrNONE, nullptr, &ltIntegral);
          if (successPropag) {
            lengthPositive = pTof.length - ltIntegral.getL();
            casctof.posFlightPr = o2::framework::pid::tof::MassToExpTime(pTof.tofExpMom, pTof.length - ltIntegral.getL(), o2::constants::physics::MassProton * o2::constants::physics::MassProton);
            casctof.posFlightPi = o2::framework::pid::tof::MassToExpTime(pTof.tofExpMom, pTof.length - ltIntegral.getL(), o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);
          }
        }
      }
      if (lengthPositive > 0.0f) {
        casctof.posDeltaTimeAsXiPi = (pTof.tofSignal - pTof.tofEvTime) - (xiFlight + lambdaFlight + casctof.posFlightPi);
        casctof.posDeltaTimeAsXiPr = (pTof.tofSignal - pTof.tofEvTime) - (xiFlight + lambdaFlight + casctof.posFlightPr);
        casctof.posDeltaTimeAsOmPi = (pTof.tofSignal - pTof.tofEvTime) - (omFlight + lambdaFlight + casctof.posFlightPi);
        casctof.posDeltaTimeAsOmPr = (pTof.tofSignal - pTof.tofEvTime) - (omFlight + lambdaFlight + casctof.posFlightPr);

        // de facto nsigma
        if (nSigmaCalibLoaded) {
          if (cascade.sign() < 0) {
            casctof.nSigmaXiLaPr = (casctof.posDeltaTimeAsXiPr - hMeanPosXiPr->Interpolate(cascade.p())) / hSigmaPosXiPr->Interpolate(cascade.p());
            casctof.nSigmaOmLaPr = (casctof.posDeltaTimeAsOmPr - hMeanPosOmPr->Interpolate(cascade.p())) / hSigmaPosOmPr->Interpolate(cascade.p());
          } else {
            casctof.nSigmaXiLaPi = (casctof.posDeltaTimeAsXiPi - hMeanPosXiPi->Interpolate(cascade.p())) / hSigmaPosXiPi->Interpolate(cascade.p());
            casctof.nSigmaOmLaPi = (casctof.posDeltaTimeAsOmPi - hMeanPosOmPi->Interpolate(cascade.p())) / hSigmaPosOmPi->Interpolate(cascade.p());
          }
        }

        // do QA histograms (calibration / QC)
        if (doQA) {
          if (passesQAcuts) {
            if (xiMinusCandidate) {
              histos.fill(HIST("h2dposDeltaTimeAsXiPr"), cascade.p(), cascade.eta(), casctof.posDeltaTimeAsXiPr);
            }
            if (xiPlusCandidate) {
              histos.fill(HIST("h2dposDeltaTimeAsXiPi"), cascade.p(), cascade.eta(), casctof.posDeltaTimeAsXiPi);
            }
            if (omegaMinusCandidate) {
              histos.fill(HIST("h2dposDeltaTimeAsOmPr"), cascade.p(), cascade.eta(), casctof.posDeltaTimeAsOmPr);
            }
            if (omegaPlusCandidate) {
              histos.fill(HIST("h2dposDeltaTimeAsOmPi"), cascade.p(), cascade.eta(), casctof.posDeltaTimeAsOmPi);
            }
          }
        }
      }
    } // end positive

    if (nTof.hasTOF && nTof.hasITS) {
      float velocityNegativePr, velocityNegativePi, lengthNegative;
      velocityNegativePr = velocityNegativePi = lengthNegative = o2::aod::v0data::kNoTOFValue;
      // method 0: legacy standalone without use of primary particle TOF
      if (calculationMethod.value == 0) {
        velocityNegativePr = velocity(negTrack.getP(), o2::constants::physics::MassProton);
        velocityNegativePi = velocity(negTrack.getP(), o2::constants::physics::MassPionCharged);
        lengthNegative = findInterceptLength(negTrack, d_bz);
        casctof.negFlightPr = lengthNegative / velocityNegativePr;
        casctof.negFlightPi = lengthNegative / velocityNegativePi;
      }
      // method 1: correct primary particle TOF information
      // length -> revise by removing travel length to primary vertex
      // expected momentum -> kept as is for now, could correct at second stage
      // use main method from TOF to calculate expected time
      if (calculationMethod.value == 1) {
        if (nTof.collisionId >= 0) {
          auto trackCollision = collisions.rawIteratorAt(nTof.collisionId);
          const o2::math_utils::Point3D<float> trackVertex{trackCollision.posX(), trackCollision.posY(), trackCollision.posZ()};
          o2::track::TrackLTIntegral ltIntegral;
          bool successPropag = o2::base::Propagator::Instance()->propagateToDCA(trackVertex, negTrack, d_bz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrNONE, nullptr, &ltIntegral);
          if (successPropag) {
            lengthNegative = nTof.length - ltIntegral.getL();
            casctof.negFlightPr = o2::framework::pid::tof::MassToExpTime(nTof.tofExpMom, nTof.length - ltIntegral.getL(), o2::constants::physics::MassProton * o2::constants::physics::MassProton);
            casctof.negFlightPi = o2::framework::pid::tof::MassToExpTime(nTof.tofExpMom, nTof.length - ltIntegral.getL(), o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);
          }
        }
      }
      if (lengthNegative > 0.0f) {
        casctof.negDeltaTimeAsXiPi = (nTof.tofSignal - nTof.tofEvTime) - (xiFlight + lambdaFlight + casctof.negFlightPi);
        casctof.negDeltaTimeAsXiPr = (nTof.tofSignal - nTof.tofEvTime) - (xiFlight + lambdaFlight + casctof.negFlightPr);
        casctof.negDeltaTimeAsOmPi = (nTof.tofSignal - nTof.tofEvTime) - (omFlight + lambdaFlight + casctof.negFlightPi);
        casctof.negDeltaTimeAsOmPr = (nTof.tofSignal - nTof.tofEvTime) - (omFlight + lambdaFlight + casctof.negFlightPr);

        // de facto nsigma
        if (nSigmaCalibLoaded) {
          if (cascade.sign() < 0) {
            casctof.nSigmaXiLaPr = (casctof.negDeltaTimeAsXiPr - hMeanPosXiPr->Interpolate(cascade.p())) / hSigmaPosXiPr->Interpolate(cascade.p());
            casctof.nSigmaOmLaPr = (casctof.negDeltaTimeAsOmPr - hMeanPosOmPr->Interpolate(cascade.p())) / hSigmaPosOmPr->Interpolate(cascade.p());
          } else {
            casctof.nSigmaXiLaPi = (casctof.negDeltaTimeAsXiPi - hMeanPosXiPi->Interpolate(cascade.p())) / hSigmaPosXiPi->Interpolate(cascade.p());
            casctof.nSigmaOmLaPi = (casctof.negDeltaTimeAsOmPi - hMeanPosOmPi->Interpolate(cascade.p())) / hSigmaPosOmPi->Interpolate(cascade.p());
          }
        }

        // do QA histograms (calibration / QC)
        if (doQA) {
          if (passesQAcuts) {
            if (xiMinusCandidate) {
              histos.fill(HIST("h2dnegDeltaTimeAsXiPi"), cascade.p(), cascade.eta(), casctof.negDeltaTimeAsXiPi);
            }
            if (xiPlusCandidate) {
              histos.fill(HIST("h2dnegDeltaTimeAsXiPr"), cascade.p(), cascade.eta(), casctof.negDeltaTimeAsXiPr);
            }
            if (omegaMinusCandidate) {
              histos.fill(HIST("h2dnegDeltaTimeAsOmPi"), cascade.p(), cascade.eta(), casctof.negDeltaTimeAsOmPi);
            }
            if (omegaPlusCandidate) {
              histos.fill(HIST("h2dnegDeltaTimeAsOmPr"), cascade.p(), cascade.eta(), casctof.negDeltaTimeAsOmPr);
            }
          }
        }
      }
    } // end negative

    if (bTof.hasTOF && bTof.hasITS) {
      float velocityBachelorKa, velocityBachelorPi, lengthBachelor;
      velocityBachelorKa = velocityBachelorPi = lengthBachelor = o2::aod::v0data::kNoTOFValue;
      // method 0: legacy standalone without use of primary particle TOF
      if (calculationMethod.value == 0) {
        velocityBachelorPi = velocity(bachTrack.getP(), o2::constants::physics::MassPionCharged);
        velocityBachelorKa = velocity(bachTrack.getP(), o2::constants::physics::MassKaonCharged);
        lengthBachelor = findInterceptLength(bachTrack, d_bz);
        casctof.bachFlightPi = lengthBachelor / velocityBachelorPi;
        casctof.bachFlightKa = lengthBachelor / velocityBachelorKa;
      }
      // method 1: correct primary particle TOF information
      // length -> revise by removing travel length to primary vertex
      // expected momentum -> kept as is for now, could correct at second stage
      // use main method from TOF to calculate expected time
      if (calculationMethod.value == 1) {
        if (bTof.collisionId >= 0) {
          auto trackCollision = collisions.rawIteratorAt(bTof.collisionId);
          const o2::math_utils::Point3D<float> trackVertex{trackCollision.posX(), trackCollision.posY(), trackCollision.posZ()};
          o2::track::TrackLTIntegral ltIntegral;
          bool successPropag = o2::base::Propagator::Instance()->propagateToDCA(trackVertex, bachTrack, d_bz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrNONE, nullptr, &ltIntegral);
          if (successPropag) {
            lengthBachelor = bTof.length - ltIntegral.getL();
            casctof.bachFlightPi = o2::framework::pid::tof::MassToExpTime(bTof.tofExpMom, bTof.length - ltIntegral.getL(), o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);
            casctof.bachFlightKa = o2::framework::pid::tof::MassToExpTime(bTof.tofExpMom, bTof.length - ltIntegral.getL(), o2::constants::physics::MassKaonCharged * o2::constants::physics::MassKaonCharged);
          }
        }
      }
      if (lengthBachelor > 0.0f) {
        casctof.bachDeltaTimeAsXiPi = (bTof.tofSignal - bTof.tofEvTime) - (xiFlight + casctof.bachFlightPi);
        casctof.bachDeltaTimeAsOmKa = (bTof.tofSignal - bTof.tofEvTime) - (omFlight + casctof.bachFlightKa);

        // de facto nsigma
        if (nSigmaCalibLoaded) {
          if (cascade.sign() < 0) {
            casctof.nSigmaXiPi = (casctof.bachDeltaTimeAsXiPi - hMeanBachXiPi->Interpolate(cascade.p())) / hSigmaBachXiPi->Interpolate(cascade.p());
            casctof.nSigmaOmKa = (casctof.bachDeltaTimeAsOmKa - hMeanBachOmKa->Interpolate(cascade.p())) / hSigmaBachOmKa->Interpolate(cascade.p());
          } else {
            casctof.nSigmaXiPi = (casctof.bachDeltaTimeAsXiPi - hMeanBachXiPi->Interpolate(cascade.p())) / hSigmaBachXiPi->Interpolate(cascade.p());
            casctof.nSigmaOmKa = (casctof.bachDeltaTimeAsOmKa - hMeanBachOmKa->Interpolate(cascade.p())) / hSigmaBachOmKa->Interpolate(cascade.p());
          }
        }

        // do QA histograms (calibration / QC)
        if (doQA) {
          if (passesQAcuts) {
            if (xiMinusCandidate) {
              histos.fill(HIST("h2dbachDeltaTimeAsXiPi"), cascade.p(), cascade.eta(), casctof.bachDeltaTimeAsXiPi);
            }
            if (xiPlusCandidate) {
              histos.fill(HIST("h2dbachDeltaTimeAsXiPi"), cascade.p(), cascade.eta(), casctof.bachDeltaTimeAsXiPi);
            }
            if (omegaMinusCandidate) {
              histos.fill(HIST("h2dbachDeltaTimeAsOmKa"), cascade.p(), cascade.eta(), casctof.bachDeltaTimeAsOmKa);
            }
            if (omegaPlusCandidate) {
              histos.fill(HIST("h2dbachDeltaTimeAsOmKa"), cascade.p(), cascade.eta(), casctof.bachDeltaTimeAsOmKa);
            }
          }
        }
      }
    } // end bachelor

    // don't forget to give feedback
    return casctof;
  }

  void processStandardData(aod::Collisions const& collisions, V0OriginalDatas const& V0s, CascOriginalDatas const& cascades, TracksWithAllExtras const& tracks, aod::BCsWithTimestamps const& /*bcs*/)
  {
    // Fire up CCDB with first collision in record. If no collisions, bypass
    if (useCustomRunNumber || collisions.size() < 1) {
      initCCDB(manualRunNumber);
    } else {
      auto collision = collisions.begin();
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc.runNumber());
    }

    //________________________________________________________________________
    // estimate event times (only necessary for original data)
    std::vector<double> collisionEventTime(collisions.size(), 0.0);
    std::vector<int> collisionNtracks(collisions.size(), 0);
    for (const auto& track : tracks) {
      if (track.hasTOF()) {
        collisionEventTime[track.collisionId()] += track.tofEvTime();
        collisionNtracks[track.collisionId()]++;
      }
    }
    for (const auto& collision : collisions) {
      if (collisionNtracks[collision.globalIndex()] > 0) {
        collisionEventTime[collision.globalIndex()] /= static_cast<double>(collisionNtracks[collision.globalIndex()]);
      } else {
        collisionEventTime[collision.globalIndex()] = -1e+6; // undefined
      }
    }

    if (calculateV0s.value) {
      for (const auto& V0 : V0s) {
        trackTofInfo pTof, nTof; // information storage

        auto pTra = V0.posTrack_as<TracksWithAllExtras>();
        auto nTra = V0.negTrack_as<TracksWithAllExtras>();

        pTof.collisionId = pTra.collisionId();
        pTof.hasITS = pTra.hasITS();
        pTof.hasTPC = pTra.hasTPC();
        pTof.hasTOF = pTra.hasTOF();
        pTof.tofExpMom = pTra.tofExpMom();
        pTof.tofEvTime = collisionEventTime[V0.collisionId()];
        pTof.tofSignal = pTra.tofSignal();
        pTof.length = pTra.length();
        pTof.tpcNSigmaPi = pTra.tpcNSigmaPi();
        pTof.tpcNSigmaPr = pTra.tpcNSigmaPr();

        nTof.collisionId = nTra.collisionId();
        nTof.hasITS = nTra.hasITS();
        nTof.hasTPC = nTra.hasTPC();
        nTof.hasTOF = nTra.hasTOF();
        nTof.tofExpMom = nTra.tofExpMom();
        nTof.tofEvTime = collisionEventTime[V0.collisionId()];
        nTof.tofSignal = nTra.tofSignal();
        nTof.length = nTra.length();
        nTof.tpcNSigmaPi = nTra.tpcNSigmaPi();
        nTof.tpcNSigmaPr = nTra.tpcNSigmaPr();

        v0TofInfo v0tof = calculateTofInfoV0(collisions, V0.collisionId(), V0, pTof, nTof);

        if (doNSigmas) {
          v0tofnsigmas(
            v0tof.nSigmaPositiveLambdaPr, v0tof.nSigmaNegativeLambdaPi,
            v0tof.nSigmaNegativeLambdaPr, v0tof.nSigmaPositiveLambdaPi,
            v0tof.nSigmaPositiveK0ShortPi, v0tof.nSigmaNegativeK0ShortPi);
        }
      }
    }

    if (calculateCascades.value) {
      for (const auto& cascade : cascades) {
        trackTofInfo pTof, nTof, bTof; // information storage

        auto pTra = cascade.posTrack_as<TracksWithAllExtras>();
        auto nTra = cascade.negTrack_as<TracksWithAllExtras>();
        auto bTra = cascade.bachelor_as<TracksWithAllExtras>();

        pTof.collisionId = pTra.collisionId();
        pTof.hasITS = pTra.hasITS();
        pTof.hasTPC = pTra.hasTPC();
        pTof.hasTOF = pTra.hasTOF();
        pTof.tofExpMom = pTra.tofExpMom();
        pTof.tofEvTime = collisionEventTime[cascade.collisionId()];
        pTof.tofSignal = pTra.tofSignal();
        pTof.length = pTra.length();
        pTof.tpcNSigmaPi = pTra.tpcNSigmaPi();
        pTof.tpcNSigmaPr = pTra.tpcNSigmaPr();

        nTof.collisionId = nTra.collisionId();
        nTof.hasITS = nTra.hasITS();
        nTof.hasTPC = nTra.hasTPC();
        nTof.hasTOF = nTra.hasTOF();
        nTof.tofExpMom = nTra.tofExpMom();
        nTof.tofEvTime = collisionEventTime[cascade.collisionId()];
        nTof.tofSignal = nTra.tofSignal();
        nTof.length = nTra.length();
        nTof.tpcNSigmaPi = nTra.tpcNSigmaPi();
        nTof.tpcNSigmaPr = nTra.tpcNSigmaPr();

        bTof.collisionId = bTra.collisionId();
        bTof.hasITS = bTra.hasITS();
        bTof.hasTPC = bTra.hasTPC();
        bTof.hasTOF = bTra.hasTOF();
        bTof.tofExpMom = bTra.tofExpMom();
        bTof.tofEvTime = collisionEventTime[cascade.collisionId()];
        bTof.tofSignal = bTra.tofSignal();
        bTof.length = bTra.length();
        bTof.tpcNSigmaPi = bTra.tpcNSigmaPi();
        bTof.tpcNSigmaKa = bTra.tpcNSigmaKa();

        cascTofInfo casctof = calculateTofInfoCascade(collisions, cascade.collisionId(), cascade, pTof, nTof, bTof);

        if (doNSigmas) {
          casctofnsigmas(
            casctof.nSigmaXiLaPi, casctof.nSigmaXiLaPr, casctof.nSigmaXiPi,
            casctof.nSigmaOmLaPi, casctof.nSigmaOmLaPr, casctof.nSigmaOmKa);
        }
      }
    }
  }

  void processDerivedData(soa::Join<aod::StraCollisions, aod::StraStamps, aod::StraEvTimes> const& collisions, V0DerivedDatas const& V0s, CascDerivedDatas const& cascades, dauTracks const& dauTrackTable, aod::DauTrackTOFPIDs const& dauTrackTOFPIDs)
  {
    // auto-determine if current or old generation of dauTrackTOFPIDs
    if (dauTrackTOFPIDs.size() == 0) {
      return;
    }
    auto firstTOFPID = dauTrackTOFPIDs.rawIteratorAt(0);
    bool isNewTOFFormat = firstTOFPID.straCollisionId() < 0 ? false : true;

    LOGF(info, "Processing derived data. Is this the new TOF info format? %i", isNewTOFFormat);

    if (!isNewTOFFormat && calculationMethod.value > 0) {
      LOGF(fatal, "Using the old derived data format with the new calculation method is not viable due to lack of needed info! Crashing.");
    }

    // Fire up CCDB with first collision in record. If no collisions, bypass
    if (useCustomRunNumber || collisions.size() < 1) {
      initCCDB(manualRunNumber);
    } else {
      auto collision = collisions.begin();
      initCCDB(collision.runNumber());
    }

    // hold indices
    std::vector<int> tofIndices(dauTrackTable.size(), -1);

    if (isNewTOFFormat) {
      // re-index
      for (const auto& dauTrackTOFPID : dauTrackTOFPIDs) {
        tofIndices[dauTrackTOFPID.dauTrackExtraId()] = dauTrackTOFPID.globalIndex();
      }
    } else {
      // they are actually joinable
      std::iota(tofIndices.begin(), tofIndices.end(), 0);
    }

    if (calculateV0s.value) {
      for (const auto& V0 : V0s) {
        trackTofInfo pTof, nTof; // information storage

        auto collision = collisions.rawIteratorAt(V0.straCollisionId());
        auto pTra = V0.posTrackExtra_as<dauTracks>();
        auto nTra = V0.negTrackExtra_as<dauTracks>();

        pTof.hasITS = pTra.hasITS();
        pTof.hasTPC = pTra.hasTPC();
        pTof.hasTOF = pTra.hasTOF();
        pTof.tpcNSigmaPi = pTra.tpcNSigmaPi();
        pTof.tpcNSigmaPr = pTra.tpcNSigmaPr();
        if (tofIndices[V0.posTrackExtraId()] >= 0 && collision.eventTime() > -1e+5) {
          auto pTofExt = dauTrackTOFPIDs.rawIteratorAt(tofIndices[V0.posTrackExtraId()]);
          pTof.collisionId = pTofExt.straCollisionId();
          pTof.tofExpMom = pTofExt.tofExpMom();
          pTof.tofEvTime = collision.eventTime();
          pTof.tofSignal = pTofExt.tofSignal();
          pTof.length = pTofExt.length();
        }

        nTof.hasITS = nTra.hasITS();
        nTof.hasTPC = nTra.hasTPC();
        nTof.hasTOF = nTra.hasTOF();
        nTof.tpcNSigmaPi = nTra.tpcNSigmaPi();
        nTof.tpcNSigmaPr = nTra.tpcNSigmaPr();
        if (tofIndices[V0.negTrackExtraId()] >= 0 && collision.eventTime() > -1e+5) {
          auto nTofExt = dauTrackTOFPIDs.rawIteratorAt(tofIndices[V0.negTrackExtraId()]);
          nTof.collisionId = nTofExt.straCollisionId();
          nTof.tofExpMom = nTofExt.tofExpMom();
          nTof.tofEvTime = collision.eventTime();
          nTof.tofSignal = nTofExt.tofSignal();
          nTof.length = nTofExt.length();
        }

        v0TofInfo v0tof = calculateTofInfoV0(collisions, V0.straCollisionId(), V0, pTof, nTof);

        if (doNSigmas) {
          v0tofnsigmas(
            v0tof.nSigmaPositiveLambdaPr, v0tof.nSigmaNegativeLambdaPi,
            v0tof.nSigmaNegativeLambdaPr, v0tof.nSigmaPositiveLambdaPi,
            v0tof.nSigmaPositiveK0ShortPi, v0tof.nSigmaNegativeK0ShortPi);
        }
      }
    }

    if (calculateCascades.value) {
      for (const auto& cascade : cascades) {
        trackTofInfo pTof, nTof, bTof; // information storage

        auto collision = collisions.rawIteratorAt(cascade.straCollisionId());
        auto pTra = cascade.posTrackExtra_as<dauTracks>();
        auto nTra = cascade.negTrackExtra_as<dauTracks>();
        auto bTra = cascade.bachTrackExtra_as<dauTracks>();

        pTof.hasITS = pTra.hasITS();
        pTof.hasTPC = pTra.hasTPC();
        pTof.hasTOF = pTra.hasTOF();
        pTof.tpcNSigmaPi = pTra.tpcNSigmaPi();
        pTof.tpcNSigmaPr = pTra.tpcNSigmaPr();
        if (tofIndices[cascade.posTrackExtraId()] >= 0 && collision.eventTime() > -1e+5) {
          auto pTofExt = dauTrackTOFPIDs.rawIteratorAt(tofIndices[cascade.posTrackExtraId()]);
          pTof.collisionId = pTofExt.straCollisionId();
          pTof.tofExpMom = pTofExt.tofExpMom();
          pTof.tofEvTime = collision.eventTime();
          pTof.tofSignal = pTofExt.tofSignal();
          pTof.length = pTofExt.length();
        }

        nTof.hasITS = nTra.hasITS();
        nTof.hasTPC = nTra.hasTPC();
        nTof.hasTOF = nTra.hasTOF();
        nTof.tpcNSigmaPi = nTra.tpcNSigmaPi();
        nTof.tpcNSigmaPr = nTra.tpcNSigmaPr();
        if (tofIndices[cascade.negTrackExtraId()] >= 0 && collision.eventTime() > -1e+5) {
          auto nTofExt = dauTrackTOFPIDs.rawIteratorAt(tofIndices[cascade.negTrackExtraId()]);
          nTof.collisionId = nTofExt.straCollisionId();
          nTof.tofExpMom = nTofExt.tofExpMom();
          nTof.tofEvTime = collision.eventTime();
          nTof.tofSignal = nTofExt.tofSignal();
          nTof.length = nTofExt.length();
        }

        bTof.hasITS = bTra.hasITS();
        bTof.hasTPC = bTra.hasTPC();
        bTof.hasTOF = bTra.hasTOF();
        bTof.tpcNSigmaPi = bTra.tpcNSigmaPi();
        bTof.tpcNSigmaKa = bTra.tpcNSigmaKa();
        if (tofIndices[cascade.bachTrackExtraId()] >= 0 && collision.eventTime() > -1e+5) {
          auto bTofExt = dauTrackTOFPIDs.rawIteratorAt(tofIndices[cascade.bachTrackExtraId()]);
          bTof.collisionId = bTofExt.straCollisionId();
          bTof.tofExpMom = bTofExt.tofExpMom();
          bTof.tofEvTime = collision.eventTime();
          bTof.tofSignal = bTofExt.tofSignal();
          bTof.length = bTofExt.length();
        }

        cascTofInfo casctof = calculateTofInfoCascade(collisions, cascade.straCollisionId(), cascade, pTof, nTof, bTof);

        if (doNSigmas) {
          casctofnsigmas(
            casctof.nSigmaXiLaPi, casctof.nSigmaXiLaPr, casctof.nSigmaXiPi,
            casctof.nSigmaOmLaPi, casctof.nSigmaOmLaPr, casctof.nSigmaOmKa);
        }
      }
    }
  }

  PROCESS_SWITCH(strangenesstofpid, processStandardData, "Process standard data", false);
  PROCESS_SWITCH(strangenesstofpid, processDerivedData, "Process derived data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenesstofpid>(cfgc)};
}
