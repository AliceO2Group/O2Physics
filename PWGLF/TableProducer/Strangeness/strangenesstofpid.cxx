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
//     Lambdakzero PID
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
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs, aod::DauTrackTOFPIDs>;
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
  Configurable<bool> correctELossInclination{"correctELossInclination", true, "factor out inclination when doing effective e-loss correction (0: no, 1: yes)"};

  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<float> tofPosition{"tofPosition", 377.934f, "TOF effective (inscribed) radius"};
  Configurable<bool> doQA{"doQA", true, "create QA histos"};
  Configurable<bool> doNSigmas{"doNSigmas", false, "calculate TOF N-sigma"};
  Configurable<bool> doQANSigma{"doQANSigma", true, "create QA of Nsigma histos"};

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

  // master p axis
  ConfigurableAxis axisP{"axisP", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};

  // for zooming in at low values only (e-loss studies and effective correction)
  ConfigurableAxis axisSmallP{"axisSmallP", {250, 0.0f, 2.5f}, "p_{T} (GeV/c)"};

  // for n-sigma calibration
  bool nSigmaCalibLoaded;
  TList* nSigmaCalibObjects;
  TH1 *hMeanPosLaPi, *hSigmaPosLaPi;
  TH1 *hMeanPosLaPr, *hSigmaPosLaPr;
  TH1 *hMeanNegLaPi, *hSigmaNegLaPi;
  TH1 *hMeanNegLaPr, *hSigmaNegLaPr;
  TH1 *hMeanPosK0Pi, *hSigmaPosK0Pi;
  TH1 *hMeanNegK0Pi, *hSigmaNegK0Pi;
  TH1 *hMeanPosXiPi, *hSigmaPosXiPi;
  TH1 *hMeanPosXiPr, *hSigmaPosXiPr;
  TH1 *hMeanNegXiPi, *hSigmaNegXiPi;
  TH1 *hMeanNegXiPr, *hSigmaNegXiPr;
  TH1 *hMeanBachXiPi, *hSigmaBachXiPi;
  TH1 *hMeanPosOmPi, *hSigmaPosOmPi;
  TH1 *hMeanPosOmPr, *hSigmaPosOmPr;
  TH1 *hMeanNegOmPi, *hSigmaNegOmPi;
  TH1 *hMeanNegOmPr, *hSigmaNegOmPr;
  TH1 *hMeanBachOmKa, *hSigmaBachOmKa;

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation

  // enum to keep track of the TOF-related properties for V0s
  enum tofEnum { kLength = 0,
                 kHasTOF,
                 kNEnums };

  /// function to calculate track length of this track up to a certain segment of a detector
  /// to be used internally in another funcrtion that calculates length until it finds the proper one
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
      float halfWidth = tofPosition * TMath::Tan(0.5f * segmentAngle);
      float x1 = TMath::Cos(theta) * (-halfWidth) + TMath::Sin(theta) * tofPosition;
      float y1 = -TMath::Sin(theta) * (-halfWidth) + TMath::Cos(theta) * tofPosition;
      float x2 = TMath::Cos(theta) * (+halfWidth) + TMath::Sin(theta) * tofPosition;
      float y2 = -TMath::Sin(theta) * (+halfWidth) + TMath::Cos(theta) * tofPosition;
      float thisLength = trackLengthToSegment(track, x1, y1, x2, y2, magneticField);
      if (thisLength < length && thisLength > 0)
        length = thisLength;
    }
    if (length > 1e+5)
      length = -100; // force negative to avoid misunderstandings
    return length;
  }

  /// O2 Propagator + TrackLTIntegral approach helpers

  /// function to calculate segmented (truncated) radius based on a certain x, y position
  float segmentedRadius(float x, float y)
  {
    float atAngle = std::atan2(y, x);
    float roundedAngle = TMath::Pi() / 9; // 18 segments = use 9 here
    float angleSegmentAxis = 0.5f * roundedAngle + roundedAngle * static_cast<float>(std::floor(atAngle / roundedAngle));
    float xSegmentAxis = TMath::Cos(angleSegmentAxis);
    float ySegmentAxis = TMath::Sin(angleSegmentAxis);
    return xSegmentAxis * x + ySegmentAxis * y; // inner product
  }

  /// function to calculate track length of this track up to a certain segmented detector
  /// \param track the input track
  /// \param time returned time (with PID given by track PID)
  void calculateTOF(o2::track::TrackPar track, float& time)
  {
    time = -1e+6;

    o2::track::TrackLTIntegral ltIntegral;

    float trackX = -100;
    static constexpr float MAX_SIN_PHI = 0.85f;
    static constexpr float MAX_STEP = 2.0f;
    static constexpr float MAX_STEP_FINAL_STAGE = 0.5f;
    static constexpr float MAX_FINAL_X = 390.0f; // maximum extra X on top of TOF X for correcting value

    bool trackOK = track.getXatLabR(tofPosition, trackX, d_bz);
    if (trackOK) {
      // propagate outwards to TOF: bulk of propagation
      o2::base::Propagator::Instance()->propagateToX(track, trackX, d_bz, MAX_SIN_PHI, MAX_STEP, o2::base::Propagator::MatCorrType::USEMatCorrLUT, &ltIntegral);

      // mark start position, define variables
      std::array<float, 3> xyz;
      track.getXYZGlo(xyz);
      float segmentedR = segmentedRadius(xyz[0], xyz[1]);
      float currentTime = ltIntegral.getTOF(track.getPID());
      if (calculationMethod.value == 2) {
        histos.fill(HIST("hTOFPosition"), xyz[0], xyz[1]); // for debugging purposes
      }

      // correct for TOF segmentation
      float trackXextra = trackX;
      bool trackOKextra = true;
      while (trackXextra < MAX_FINAL_X) {
        // propagate one step further
        trackXextra += MAX_STEP_FINAL_STAGE;
        trackOKextra = o2::base::Propagator::Instance()->propagateToX(track, trackXextra, d_bz, MAX_SIN_PHI, MAX_STEP, o2::base::Propagator::MatCorrType::USEMatCorrLUT, &ltIntegral);
        if (!trackOKextra) {
          time = -1e+6;
          return; // propagation failed, skip, won't look reasonable
        }

        // re-evaluate - did we cross? if yes break
        float previousX = xyz[0], previousY = xyz[1];
        track.getXYZGlo(xyz);
        if (segmentedRadius(xyz[0], xyz[1]) > tofPosition) {
          // crossed boundary -> do proportional scaling with how much we actually crossed the boundary
          float segmentedRFinal = segmentedRadius(xyz[0], xyz[1]);
          float timeFinal = ltIntegral.getTOF(track.getPID());
          float fraction = (tofPosition - segmentedR) / (segmentedRFinal - segmentedR + 1e-6); // proportional fraction
          time = currentTime + (timeFinal - currentTime) * fraction;
          if (calculationMethod.value == 2) {
            histos.fill(HIST("hTOFPositionFinal"), previousX + fraction * (xyz[0] - previousX), previousY + fraction * (xyz[1] - previousY)); // for debugging purposes
          }
          return; // get out of the entire function and return (don't just break)
        }

        // prepare for next step by setting current position and desired variables
        segmentedR = segmentedRadius(xyz[0], xyz[1]);
        currentTime = ltIntegral.getTOF(track.getPID());
      }
    }
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

    if (calculationMethod.value == 2) {
      //_____________________________________________________________________
      // special mode in which comparison histograms are required

      // base ArcDebug: comparison between times of arrival in different methods
      histos.add("hArcDebug", "hArcDebug", kTH2F, {axisTime, axisTime});

      // Position of TrackLTIntegral method: intermediate (getXatLabR) and final (reach segmented detector)
      histos.add("hTOFPosition", "hTOFPosition", kTH2F, {axisPosition, axisPosition});
      histos.add("hTOFPositionFinal", "hTOFPositionFinal", kTH2F, {axisPosition, axisPosition});

      // Delta-times of each method for the various species
      histos.add("hDeltaTimeMethodsVsP_posLaPr", "hDeltaTimeMethodsVsP_posLaPr", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_posLaPi", "hDeltaTimeMethodsVsP_posLaPi", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_posK0Pi", "hDeltaTimeMethodsVsP_posK0Pi", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_negLaPr", "hDeltaTimeMethodsVsP_negLaPr", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_negLaPi", "hDeltaTimeMethodsVsP_negLaPi", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_negK0Pi", "hDeltaTimeMethodsVsP_negK0Pi", kTH3F, {axisSmallP, axisEta, axisDeltaTime});

      histos.add("hDeltaTimeMethodsVsP_posXiPi", "hDeltaTimeMethodsVsP_posXiPi", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_posXiPr", "hDeltaTimeMethodsVsP_posXiPr", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_negXiPi", "hDeltaTimeMethodsVsP_negXiPi", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_negXiPr", "hDeltaTimeMethodsVsP_negXiPr", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_bachXiPi", "hDeltaTimeMethodsVsP_bachXiPi", kTH3F, {axisSmallP, axisEta, axisDeltaTime});

      histos.add("hDeltaTimeMethodsVsP_posOmPi", "hDeltaTimeMethodsVsP_posOmPi", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_posOmPr", "hDeltaTimeMethodsVsP_posOmPr", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_negOmPi", "hDeltaTimeMethodsVsP_negOmPi", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_negOmPr", "hDeltaTimeMethodsVsP_negOmPr", kTH3F, {axisSmallP, axisEta, axisDeltaTime});
      histos.add("hDeltaTimeMethodsVsP_bachOmKa", "hDeltaTimeMethodsVsP_bachOmKa", kTH3F, {axisSmallP, axisEta, axisDeltaTime});

      histos.add("hRatioTimeMethodsVsP_posLaPr", "hRatioTimeMethodsVsP_posLaPr", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_posLaPi", "hRatioTimeMethodsVsP_posLaPi", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_posK0Pi", "hRatioTimeMethodsVsP_posK0Pi", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_negLaPr", "hRatioTimeMethodsVsP_negLaPr", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_negLaPi", "hRatioTimeMethodsVsP_negLaPi", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_negK0Pi", "hRatioTimeMethodsVsP_negK0Pi", kTH3F, {axisSmallP, axisEta, axisRatioMethods});

      histos.add("hRatioTimeMethodsVsP_posXiPi", "hRatioTimeMethodsVsP_posXiPi", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_posXiPr", "hRatioTimeMethodsVsP_posXiPr", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_negXiPi", "hRatioTimeMethodsVsP_negXiPi", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_negXiPr", "hRatioTimeMethodsVsP_negXiPr", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_bachXiPi", "hRatioTimeMethodsVsP_bachXiPi", kTH3F, {axisSmallP, axisEta, axisRatioMethods});

      histos.add("hRatioTimeMethodsVsP_posOmPi", "hRatioTimeMethodsVsP_posOmPi", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_posOmPr", "hRatioTimeMethodsVsP_posOmPr", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_negOmPi", "hRatioTimeMethodsVsP_negOmPi", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_negOmPr", "hRatioTimeMethodsVsP_negOmPr", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
      histos.add("hRatioTimeMethodsVsP_bachOmKa", "hRatioTimeMethodsVsP_bachOmKa", kTH3F, {axisSmallP, axisEta, axisRatioMethods});
    }

    // list memory consumption at start if running in modes with more output
    if (calculationMethod.value == 2 || doQA) {
      histos.print();
    }
  }

  void initCCDB(int runNumber)
  {
    if (mRunNumber == runNumber) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
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

    if (calculationMethod.value > 0 && !lut) {
      // setMatLUT only after magfield has been initalized
      // (setMatLUT has implicit and problematic init field call if not)
      LOG(info) << "Loading full (all-radius) material look-up table for run number: " << runNumber;
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->getForRun<o2::base::MatLayerCylSet>(ccdbConfigurations.lutPath, runNumber));
      o2::base::Propagator::Instance()->setMatLUT(lut);
      LOG(info) << "Material look-up table loaded!";
    }
    mRunNumber = runNumber;
  }

  float velocity(float lMomentum, float lMass)
  {
    // Momentum p and mass m -> returns speed in centimeters per picosecond
    // Useful for TOF calculations
    float lA = (lMomentum / lMass) * (lMomentum / lMass);
    return 0.0299792458 * TMath::Sqrt(lA / (1 + lA));
  }

  // templatized process function for symmetric operation in derived and original AO2D
  template <class TCollision, typename TV0, typename TTrack>
  void processV0Candidate(TCollision const& collision, TV0 const& v0, TTrack const& pTra, TTrack const& nTra, int v0pdg)
  {
    // time of V0 segment
    float lengthV0 = std::hypot(v0.x() - collision.getX(), v0.y() - collision.getY(), v0.z() - collision.getZ());
    float velocityK0Short = velocity(v0.p(), o2::constants::physics::MassKaonNeutral);
    float velocityLambda = velocity(v0.p(), o2::constants::physics::MassLambda);
    float timeK0Short = lengthV0 / velocityK0Short; // in picoseconds
    float timeLambda = lengthV0 / velocityLambda;   // in picoseconds

    // initialize from V0 position and momenta
    o2::track::TrackPar posTrack = o2::track::TrackPar({v0.x(), v0.y(), v0.z()}, {v0.pxpos(), v0.pypos(), v0.pzpos()}, +1);
    o2::track::TrackPar negTrack = o2::track::TrackPar({v0.x(), v0.y(), v0.z()}, {v0.pxneg(), v0.pyneg(), v0.pzneg()}, -1);

    // at minimum
    float positiveP = std::hypot(v0.pxpos(), v0.pypos(), v0.pzpos());
    float negativeP = std::hypot(v0.pxneg(), v0.pyneg(), v0.pzneg());

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

    float timePositivePr = o2::aod::v0data::kNoTOFValue;
    float timePositivePi = o2::aod::v0data::kNoTOFValue;
    float timeNegativePr = o2::aod::v0data::kNoTOFValue;
    float timeNegativePi = o2::aod::v0data::kNoTOFValue;

    float timePositivePr_Method0 = o2::aod::v0data::kNoTOFValue;
    float timePositivePi_Method0 = o2::aod::v0data::kNoTOFValue;
    float timeNegativePr_Method0 = o2::aod::v0data::kNoTOFValue;
    float timeNegativePi_Method0 = o2::aod::v0data::kNoTOFValue;

    float timePositivePr_Method1 = o2::aod::v0data::kNoTOFValue;
    float timePositivePi_Method1 = o2::aod::v0data::kNoTOFValue;
    float timeNegativePr_Method1 = o2::aod::v0data::kNoTOFValue;
    float timeNegativePi_Method1 = o2::aod::v0data::kNoTOFValue;

    if (calculationMethod.value == 0 || calculationMethod.value == 2) {
      float velocityPositivePr = velocity(posTrack.getP(), o2::constants::physics::MassProton);
      float velocityPositivePi = velocity(posTrack.getP(), o2::constants::physics::MassPionCharged);
      float velocityNegativePr = velocity(negTrack.getP(), o2::constants::physics::MassProton);
      float velocityNegativePi = velocity(negTrack.getP(), o2::constants::physics::MassPionCharged);

      float lengthPositive = findInterceptLength(posTrack, d_bz); // FIXME: tofPosition ok? adjust?
      float lengthNegative = findInterceptLength(negTrack, d_bz); // FIXME: tofPosition ok? adjust?

      if (lengthPositive > 0) {
        timePositivePr_Method0 = lengthPositive / velocityPositivePr;
        timePositivePi_Method0 = lengthPositive / velocityPositivePi;
      }
      if (lengthNegative > 0) {
        timeNegativePr_Method0 = lengthNegative / velocityNegativePr;
        timeNegativePi_Method0 = lengthNegative / velocityNegativePi;
      }
    }

    if (calculationMethod.value > 0) {
      // method to calculate the time and length via Propagator TrackLTIntegral
      if (pTra.hasTOF()) { // calculate if signal present, otherwise skip
        o2::track::TrackPar posTrackAsProton(posTrack);
        posTrackAsProton.setPID(o2::track::PID::Proton);
        calculateTOF(posTrackAsProton, timePositivePr_Method1);

        o2::track::TrackPar posTrackAsPion(posTrack);
        posTrackAsPion.setPID(o2::track::PID::Pion);
        calculateTOF(posTrackAsPion, timePositivePi_Method1);
      }
      if (nTra.hasTOF()) { // calculate if signal present, otherwise skip
        o2::track::TrackPar negTrackAsProton(negTrack);
        negTrackAsProton.setPID(o2::track::PID::Proton);
        calculateTOF(negTrackAsProton, timeNegativePr_Method1);

        o2::track::TrackPar negTrackAsPion(negTrack);
        negTrackAsPion.setPID(o2::track::PID::Pion);
        calculateTOF(negTrackAsPion, timeNegativePi_Method1);
      }
    }

    // assign values to be used in main calculation
    if (calculationMethod.value == 0) {
      timePositivePr = timePositivePr_Method0;
      timePositivePi = timePositivePi_Method0;
      timeNegativePr = timeNegativePr_Method0;
      timeNegativePi = timeNegativePi_Method0;
    } else {
      timePositivePr = timePositivePr_Method1;
      timePositivePi = timePositivePi_Method1;
      timeNegativePr = timeNegativePr_Method1;
      timeNegativePi = timeNegativePi_Method1;
    }

    if (pTra.hasTOF() && timePositivePr > 0) {
      deltaTimePositiveLambdaPr = (pTra.tofSignal() - pTra.tofEvTime()) - (timeLambda + timePositivePr);
      deltaTimePositiveLambdaPi = (pTra.tofSignal() - pTra.tofEvTime()) - (timeLambda + timePositivePi);
      deltaTimePositiveK0ShortPi = (pTra.tofSignal() - pTra.tofEvTime()) - (timeK0Short + timePositivePi);
    }
    if (nTra.hasTOF() && timeNegativePr > 0) {
      deltaTimeNegativeLambdaPr = (nTra.tofSignal() - nTra.tofEvTime()) - (timeLambda + timeNegativePr);
      deltaTimeNegativeLambdaPi = (nTra.tofSignal() - nTra.tofEvTime()) - (timeLambda + timeNegativePi);
      deltaTimeNegativeK0ShortPi = (nTra.tofSignal() - nTra.tofEvTime()) - (timeK0Short + timeNegativePi);
    }

    if (doQA) {
      // calculate and pack properties for QA purposes
      int posProperties = 0;
      if (timePositivePr > 0)
        posProperties = posProperties | (static_cast<int>(1) << kLength);
      if (pTra.hasTOF())
        posProperties = posProperties | (static_cast<int>(1) << kHasTOF);
      int negProperties = 0;
      if (timeNegativePr > 0)
        negProperties = negProperties | (static_cast<int>(1) << kLength);
      if (nTra.hasTOF())
        negProperties = negProperties | (static_cast<int>(1) << kHasTOF);

      histos.fill(HIST("h2dPositiveTOFProperties"), v0.p(), posProperties);
      histos.fill(HIST("h2dNegativeTOFProperties"), v0.p(), negProperties);
    }

    float deltaDecayTimeLambda = -10e+4;
    float deltaDecayTimeAntiLambda = -10e+4;
    float deltaDecayTimeK0Short = -10e+4;
    if (nTra.hasTOF() && pTra.hasTOF() > 0 && timePositivePr > 0 && timeNegativePr > 0) { // does not depend on event time
      deltaDecayTimeLambda = (pTra.tofSignal() - timePositivePr) - (nTra.tofSignal() - timeNegativePi);
      deltaDecayTimeAntiLambda = (pTra.tofSignal() - timePositivePi) - (nTra.tofSignal() - timeNegativePr);
      deltaDecayTimeK0Short = (pTra.tofSignal() - timePositivePi) - (nTra.tofSignal() - timeNegativePi);
    }

    // calculate betas

    float evTimeMean = 0.5f * (pTra.tofEvTime() + nTra.tofEvTime());
    float decayTimeLambda = 0.5f * ((pTra.tofSignal() - timePositivePr) + (nTra.tofSignal() - timeNegativePi)) - evTimeMean;
    float decayTimeAntiLambda = 0.5f * ((pTra.tofSignal() - timePositivePi) + (nTra.tofSignal() - timeNegativePr)) - evTimeMean;
    float decayTimeK0Short = 0.5f * ((pTra.tofSignal() - timePositivePi) + (nTra.tofSignal() - timeNegativePi)) - evTimeMean;

    float betaLambda = o2::aod::cascdata::kNoTOFValue;
    float betaAntiLambda = o2::aod::cascdata::kNoTOFValue;
    float betaK0Short = o2::aod::cascdata::kNoTOFValue;

    if (nTra.hasTOF() && pTra.hasTOF()) {
      betaLambda = (lengthV0 / decayTimeLambda) / 0.0299792458;
      betaAntiLambda = (lengthV0 / decayTimeAntiLambda) / 0.0299792458;
      betaK0Short = (lengthV0 / decayTimeK0Short) / 0.0299792458;
    }

    v0tofpid(deltaTimePositiveLambdaPi, deltaTimePositiveLambdaPr,
             deltaTimeNegativeLambdaPi, deltaTimeNegativeLambdaPr,
             deltaTimePositiveK0ShortPi, deltaTimeNegativeK0ShortPi,
             deltaDecayTimeLambda, deltaDecayTimeAntiLambda, deltaDecayTimeK0Short);
    v0tofbeta(betaLambda, betaAntiLambda, betaK0Short);
    v0tofdebugs(timeLambda, timeK0Short, timePositivePr, timePositivePi, timeNegativePr, timeNegativePi);

    // do Nsigmas if requested
    if (doNSigmas && nSigmaCalibLoaded) {
      // sweep through all viable hypotheses and produce N-sigma

      if (deltaTimePositiveLambdaPi > -1e+5)
        nSigmaPositiveLambdaPi = (deltaTimePositiveLambdaPi - hMeanPosLaPi->Interpolate(v0.p())) / hSigmaPosLaPi->Interpolate(v0.p());
      if (deltaTimePositiveLambdaPr > -1e+5)
        nSigmaPositiveLambdaPr = (deltaTimePositiveLambdaPr - hMeanPosLaPr->Interpolate(v0.p())) / hSigmaPosLaPr->Interpolate(v0.p());
      if (deltaTimeNegativeLambdaPi > -1e+5)
        nSigmaNegativeLambdaPi = (deltaTimeNegativeLambdaPi - hMeanNegLaPi->Interpolate(v0.p())) / hSigmaNegLaPi->Interpolate(v0.p());
      if (deltaTimeNegativeLambdaPr > -1e+5)
        nSigmaNegativeLambdaPr = (deltaTimeNegativeLambdaPr - hMeanNegLaPr->Interpolate(v0.p())) / hSigmaNegLaPr->Interpolate(v0.p());
      if (deltaTimePositiveK0ShortPi > -1e+5)
        nSigmaPositiveK0ShortPi = (deltaTimePositiveK0ShortPi - hMeanPosK0Pi->Interpolate(v0.p())) / hSigmaPosK0Pi->Interpolate(v0.p());
      if (deltaTimeNegativeK0ShortPi > -1e+5)
        nSigmaNegativeK0ShortPi = (deltaTimeNegativeK0ShortPi - hMeanNegK0Pi->Interpolate(v0.p())) / hSigmaNegK0Pi->Interpolate(v0.p());

      v0tofnsigmas(
        nSigmaPositiveLambdaPr, nSigmaNegativeLambdaPi,
        nSigmaNegativeLambdaPr, nSigmaPositiveLambdaPi,
        nSigmaPositiveK0ShortPi, nSigmaNegativeK0ShortPi);
    }

    if (doQA) {
      // length factor due to eta (to offset e-loss)
      float positiveCosine = 1.0f / sqrt(1.0f + posTrack.getTgl() * posTrack.getTgl());
      float negativeCosine = 1.0f / sqrt(1.0f + negTrack.getTgl() * negTrack.getTgl());
      if (correctELossInclination.value == false) {
        negativeCosine = positiveCosine = 1.0f;
      }

      if (pTra.hasTOF()) {
        if (v0.v0cosPA() > v0Group.qaCosPA && v0.dcaV0daughters() < v0Group.qaDCADau) {
          if (std::abs(v0.mLambda() - 1.115683) < v0Group.qaMassWindow && fabs(pTra.tpcNSigmaPr()) < v0Group.qaTPCNSigma && fabs(nTra.tpcNSigmaPi()) < v0Group.qaTPCNSigma && ((v0pdg == 0) || (v0pdg == 3122))) {
            histos.fill(HIST("h2dDeltaTimePositiveLambdaPr"), v0.p(), v0.eta(), deltaTimePositiveLambdaPr);
            if (calculationMethod.value == 2 && std::abs(timePositivePr_Method0 - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon && std::abs(timePositivePr_Method1 - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
              histos.fill(HIST("hDeltaTimeMethodsVsP_posLaPr"), positiveP, v0.positiveeta(), (timePositivePr_Method0 - timePositivePr_Method1) * positiveCosine);
              histos.fill(HIST("hRatioTimeMethodsVsP_posLaPr"), positiveP, v0.positiveeta(), (timePositivePr_Method1 / timePositivePr_Method0) * positiveCosine);
            }
            if (doQANSigma)
              histos.fill(HIST("h2dNSigmaPositiveLambdaPr"), v0.p(), nSigmaPositiveLambdaPr);
          }
          if (std::abs(v0.mAntiLambda() - 1.115683) < v0Group.qaMassWindow && fabs(pTra.tpcNSigmaPi()) < v0Group.qaTPCNSigma && fabs(nTra.tpcNSigmaPr()) < v0Group.qaTPCNSigma && ((v0pdg == 0) || (v0pdg == -3122))) {
            histos.fill(HIST("h2dDeltaTimePositiveLambdaPi"), v0.p(), v0.eta(), deltaTimePositiveLambdaPi);
            if (calculationMethod.value == 2 && std::abs(timePositivePi_Method0 - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon && std::abs(timePositivePi_Method1 - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
              histos.fill(HIST("hDeltaTimeMethodsVsP_posLaPi"), positiveP, v0.positiveeta(), (timePositivePi_Method0 - timePositivePi_Method1) * positiveCosine);
              histos.fill(HIST("hRatioTimeMethodsVsP_posLaPi"), positiveP, v0.positiveeta(), (timePositivePi_Method1 / timePositivePi_Method0) * positiveCosine);
            }
            if (doQANSigma)
              histos.fill(HIST("h2dNSigmaPositiveLambdaPi"), v0.p(), nSigmaPositiveLambdaPi);
          }
          if (std::abs(v0.mK0Short() - 0.497) < v0Group.qaMassWindow && fabs(pTra.tpcNSigmaPi()) < v0Group.qaTPCNSigma && fabs(nTra.tpcNSigmaPi()) < v0Group.qaTPCNSigma && ((v0pdg == 0) || (v0pdg == 310))) {
            histos.fill(HIST("h2dDeltaTimePositiveK0ShortPi"), v0.p(), v0.eta(), deltaTimePositiveK0ShortPi);
            if (calculationMethod.value == 2 && std::abs(timePositivePi_Method0 - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon && std::abs(timePositivePi_Method1 - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
              histos.fill(HIST("hDeltaTimeMethodsVsP_posK0Pi"), positiveP, v0.positiveeta(), (timePositivePi_Method0 - timePositivePi_Method1) * positiveCosine);
              histos.fill(HIST("hRatioTimeMethodsVsP_posK0Pi"), positiveP, v0.positiveeta(), (timePositivePi_Method1 / timePositivePi_Method0) * positiveCosine);
            }
            if (doQANSigma)
              histos.fill(HIST("h2dNSigmaPositiveK0ShortPi"), v0.p(), nSigmaPositiveK0ShortPi);
          }
        }
      }

      if (nTra.hasTOF()) {
        if (v0.v0cosPA() > v0Group.qaCosPA && v0.dcaV0daughters() < v0Group.qaDCADau) {
          if (std::abs(v0.mLambda() - 1.115683) < v0Group.qaMassWindow && fabs(pTra.tpcNSigmaPr()) < v0Group.qaTPCNSigma && fabs(nTra.tpcNSigmaPi()) < v0Group.qaTPCNSigma && ((v0pdg == 0) || (v0pdg == 3122))) {
            histos.fill(HIST("h2dDeltaTimeNegativeLambdaPi"), v0.p(), v0.eta(), deltaTimeNegativeLambdaPi);
            if (calculationMethod.value == 2 && std::abs(timeNegativePi_Method0 - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon && std::abs(timeNegativePi_Method1 - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
              histos.fill(HIST("hDeltaTimeMethodsVsP_negLaPi"), negativeP, v0.negativeeta(), (timeNegativePi_Method0 - timeNegativePi_Method1) * negativeCosine);
              histos.fill(HIST("hRatioTimeMethodsVsP_negLaPi"), negativeP, v0.negativeeta(), (timeNegativePi_Method1 / timeNegativePi_Method0) * negativeCosine);
            }
            if (doQANSigma)
              histos.fill(HIST("h2dNSigmaNegativeLambdaPi"), v0.p(), nSigmaNegativeLambdaPi);
          }
          if (std::abs(v0.mAntiLambda() - 1.115683) < v0Group.qaMassWindow && fabs(pTra.tpcNSigmaPi()) < v0Group.qaTPCNSigma && fabs(nTra.tpcNSigmaPr()) < v0Group.qaTPCNSigma && ((v0pdg == 0) || (v0pdg == -3122))) {
            histos.fill(HIST("h2dDeltaTimeNegativeLambdaPr"), v0.p(), v0.eta(), deltaTimeNegativeLambdaPr);
            if (calculationMethod.value == 2 && std::abs(timeNegativePr_Method0 - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon && std::abs(timeNegativePr_Method1 - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
              histos.fill(HIST("hDeltaTimeMethodsVsP_negLaPr"), negativeP, v0.negativeeta(), (timeNegativePr_Method0 - timeNegativePr_Method1) * negativeCosine);
              histos.fill(HIST("hRatioTimeMethodsVsP_negLaPr"), negativeP, v0.negativeeta(), (timeNegativePr_Method1 / timeNegativePr_Method0) * negativeCosine);
            }
            if (doQANSigma)
              histos.fill(HIST("h2dNSigmaNegativeLambdaPr"), v0.p(), nSigmaNegativeLambdaPr);
          }
          if (std::abs(v0.mK0Short() - 0.497) < v0Group.qaMassWindow && fabs(pTra.tpcNSigmaPi()) < v0Group.qaTPCNSigma && fabs(nTra.tpcNSigmaPi()) < v0Group.qaTPCNSigma && ((v0pdg == 0) || (v0pdg == 310))) {
            histos.fill(HIST("h2dDeltaTimeNegativeK0ShortPi"), v0.p(), v0.eta(), deltaTimeNegativeK0ShortPi);
            if (calculationMethod.value == 2 && std::abs(timeNegativePi_Method0 - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon && std::abs(timeNegativePi_Method1 - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
              histos.fill(HIST("hDeltaTimeMethodsVsP_negK0Pi"), negativeP, v0.negativeeta(), (timeNegativePi_Method0 - timeNegativePi_Method1) * negativeCosine);
              histos.fill(HIST("hRatioTimeMethodsVsP_negK0Pi"), negativeP, v0.negativeeta(), (timeNegativePi_Method1 / timeNegativePi_Method0) * negativeCosine);
            }
            if (doQANSigma)
              histos.fill(HIST("h2dNSigmaNegativeK0ShortPi"), v0.p(), nSigmaNegativeK0ShortPi);
          }
        }
      }
      // delta lambda decay time
      histos.fill(HIST("h2dLambdaDeltaDecayTime"), v0.p(), deltaDecayTimeLambda);
    }
  }

  template <class TCollision, typename TCascade, typename TTrack>
  void processCascadeCandidate(TCollision const& collision, TCascade const& cascade, TTrack const& pTra, TTrack const& nTra, TTrack const& bTra, int cascpdg)
  {
    // initialize from positions and momenta as needed
    o2::track::TrackPar posTrack = o2::track::TrackPar({cascade.xlambda(), cascade.ylambda(), cascade.zlambda()}, {cascade.pxpos(), cascade.pypos(), cascade.pzpos()}, +1);
    o2::track::TrackPar negTrack = o2::track::TrackPar({cascade.xlambda(), cascade.ylambda(), cascade.zlambda()}, {cascade.pxneg(), cascade.pyneg(), cascade.pzneg()}, -1);
    o2::track::TrackPar bachTrack = o2::track::TrackPar({cascade.x(), cascade.y(), cascade.z()}, {cascade.pxbach(), cascade.pybach(), cascade.pzbach()}, cascade.sign());
    o2::track::TrackPar cascTrack = o2::track::TrackPar({cascade.x(), cascade.y(), cascade.z()}, {cascade.px(), cascade.py(), cascade.pz()}, cascade.sign());

    float positiveP = std::hypot(cascade.pxpos(), cascade.pypos(), cascade.pzpos());
    float negativeP = std::hypot(cascade.pxneg(), cascade.pyneg(), cascade.pzneg());
    float bachelorP = std::hypot(cascade.pxbach(), cascade.pybach(), cascade.pzbach());

    // start calculation: calculate velocities
    float velocityXi = velocity(cascTrack.getP(), o2::constants::physics::MassXiMinus);
    float velocityOm = velocity(cascTrack.getP(), o2::constants::physics::MassOmegaMinus);
    float velocityLa = velocity(std::hypot(cascade.pxlambda(), cascade.pylambda(), cascade.pzlambda()), o2::constants::physics::MassLambda);

    // calculate mother lengths
    float lengthV0 = std::hypot(cascade.xlambda() - cascade.x(), cascade.ylambda() - cascade.y(), cascade.zlambda() - cascade.z());
    float lengthCascade = o2::aod::cascdata::kNoTOFValue;
    ;
    const o2::math_utils::Point3D<float> collVtx{collision.getX(), collision.getY(), collision.getZ()};
    bool successPropag = o2::base::Propagator::Instance()->propagateToDCA(collVtx, cascTrack, d_bz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrNONE);
    float d = -1.0f;
    float linearToPV = std::hypot(cascade.x() - collision.getX(), cascade.y() - collision.getY(), cascade.z() - collision.getZ());
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

    // lambda, xi and omega flight time is always defined
    float lambdaFlight = lengthV0 / velocityLa;
    float xiFlight = lengthCascade / velocityXi;
    float omFlight = lengthCascade / velocityOm;
    float posFlightPi = o2::aod::cascdata::kNoTOFValue;
    float posFlightPr = o2::aod::cascdata::kNoTOFValue;
    float negFlightPi = o2::aod::cascdata::kNoTOFValue;
    float negFlightPr = o2::aod::cascdata::kNoTOFValue;
    float bachFlightPi = o2::aod::cascdata::kNoTOFValue;
    float bachFlightKa = o2::aod::cascdata::kNoTOFValue;

    float posFlightPi_Method0 = o2::aod::cascdata::kNoTOFValue;
    float posFlightPr_Method0 = o2::aod::cascdata::kNoTOFValue;
    float negFlightPi_Method0 = o2::aod::cascdata::kNoTOFValue;
    float negFlightPr_Method0 = o2::aod::cascdata::kNoTOFValue;
    float bachFlightPi_Method0 = o2::aod::cascdata::kNoTOFValue;
    float bachFlightKa_Method0 = o2::aod::cascdata::kNoTOFValue;

    float posFlightPi_Method1 = o2::aod::cascdata::kNoTOFValue;
    float posFlightPr_Method1 = o2::aod::cascdata::kNoTOFValue;
    float negFlightPi_Method1 = o2::aod::cascdata::kNoTOFValue;
    float negFlightPr_Method1 = o2::aod::cascdata::kNoTOFValue;
    float bachFlightPi_Method1 = o2::aod::cascdata::kNoTOFValue;
    float bachFlightKa_Method1 = o2::aod::cascdata::kNoTOFValue;

    // actual time-of-flight of daughter calculation
    if (calculationMethod.value == 0 || calculationMethod.value == 2) {
      float velocityPositivePr = velocity(posTrack.getP(), o2::constants::physics::MassProton);
      float velocityPositivePi = velocity(posTrack.getP(), o2::constants::physics::MassPionCharged);
      float velocityNegativePr = velocity(negTrack.getP(), o2::constants::physics::MassProton);
      float velocityNegativePi = velocity(negTrack.getP(), o2::constants::physics::MassPionCharged);
      float velocityBachelorPi = velocity(bachTrack.getP(), o2::constants::physics::MassPionCharged);
      float velocityBachelorKa = velocity(bachTrack.getP(), o2::constants::physics::MassKaonCharged);

      float lengthPositive = findInterceptLength(posTrack, d_bz);  // FIXME: tofPosition ok? adjust?
      float lengthNegative = findInterceptLength(negTrack, d_bz);  // FIXME: tofPosition ok? adjust?
      float lengthBachelor = findInterceptLength(bachTrack, d_bz); // FIXME: tofPosition ok? adjust?

      if (lengthPositive > 0) {
        posFlightPi_Method0 = lengthPositive / velocityPositivePi;
        posFlightPr_Method0 = lengthPositive / velocityPositivePr;
      }
      if (lengthNegative > 0) {
        negFlightPi_Method0 = lengthNegative / velocityNegativePi;
        negFlightPr_Method0 = lengthNegative / velocityNegativePr;
      }
      if (lengthBachelor > 0) {
        bachFlightPi_Method0 = lengthBachelor / velocityBachelorPi;
        bachFlightKa_Method0 = lengthBachelor / velocityBachelorKa;
      }
    }

    if (calculationMethod.value > 0) {
      if (pTra.hasTOF()) { // calculate if signal present, otherwise skip
        o2::track::TrackPar posTrackAsProton(posTrack);
        posTrackAsProton.setPID(o2::track::PID::Proton);
        calculateTOF(posTrackAsProton, posFlightPr_Method1);

        o2::track::TrackPar posTrackAsPion(posTrack);
        posTrackAsPion.setPID(o2::track::PID::Pion);
        calculateTOF(posTrackAsPion, posFlightPi_Method1);
      }
      if (nTra.hasTOF()) { // calculate if signal present, otherwise skip
        o2::track::TrackPar negTrackAsProton(negTrack);
        negTrackAsProton.setPID(o2::track::PID::Proton);
        calculateTOF(negTrackAsProton, negFlightPr_Method1);

        o2::track::TrackPar negTrackAsPion(negTrack);
        negTrackAsPion.setPID(o2::track::PID::Pion);
        calculateTOF(negTrackAsPion, negFlightPi_Method1);
      }
      if (bTra.hasTOF()) { // calculate if signal present, otherwise skip
        o2::track::TrackPar bachTrackAsPion(bachTrack);
        bachTrackAsPion.setPID(o2::track::PID::Pion);
        calculateTOF(bachTrackAsPion, bachFlightPi_Method1);

        o2::track::TrackPar bachTrackAsKaon(bachTrack);
        bachTrackAsKaon.setPID(o2::track::PID::Kaon);
        calculateTOF(bachTrackAsKaon, bachFlightKa_Method1);
      }
    }

    // assign values to be used in main calculation
    if (calculationMethod.value == 0) {
      posFlightPi = posFlightPi_Method0;
      posFlightPr = posFlightPr_Method0;
      negFlightPi = negFlightPi_Method0;
      negFlightPr = negFlightPr_Method0;
      bachFlightPi = bachFlightPi_Method0;
      bachFlightKa = bachFlightKa_Method0;
    } else {
      posFlightPi = posFlightPi_Method1;
      posFlightPr = posFlightPr_Method1;
      negFlightPi = negFlightPi_Method1;
      negFlightPr = negFlightPr_Method1;
      bachFlightPi = bachFlightPi_Method1;
      bachFlightKa = bachFlightKa_Method1;
    }

    // initialize delta-times (actual PID variables)
    float posDeltaTimeAsXiPi = o2::aod::cascdata::kNoTOFValue, posDeltaTimeAsXiPr = o2::aod::cascdata::kNoTOFValue;
    float negDeltaTimeAsXiPi = o2::aod::cascdata::kNoTOFValue, negDeltaTimeAsXiPr = o2::aod::cascdata::kNoTOFValue;
    float bachDeltaTimeAsXiPi = o2::aod::cascdata::kNoTOFValue;
    float posDeltaTimeAsOmPi = o2::aod::cascdata::kNoTOFValue, posDeltaTimeAsOmPr = o2::aod::cascdata::kNoTOFValue;
    float negDeltaTimeAsOmPi = o2::aod::cascdata::kNoTOFValue, negDeltaTimeAsOmPr = o2::aod::cascdata::kNoTOFValue;
    float bachDeltaTimeAsOmKa = o2::aod::cascdata::kNoTOFValue;

    if (pTra.hasTOF()) {
      posDeltaTimeAsXiPi = (pTra.tofSignal() - pTra.tofEvTime()) - (xiFlight + lambdaFlight + posFlightPi);
      posDeltaTimeAsXiPr = (pTra.tofSignal() - pTra.tofEvTime()) - (xiFlight + lambdaFlight + posFlightPr);
      posDeltaTimeAsOmPi = (pTra.tofSignal() - pTra.tofEvTime()) - (omFlight + lambdaFlight + posFlightPi);
      posDeltaTimeAsOmPr = (pTra.tofSignal() - pTra.tofEvTime()) - (omFlight + lambdaFlight + posFlightPr);
    }
    if (nTra.hasTOF()) {
      negDeltaTimeAsXiPi = (nTra.tofSignal() - nTra.tofEvTime()) - (xiFlight + lambdaFlight + negFlightPi);
      negDeltaTimeAsXiPr = (nTra.tofSignal() - nTra.tofEvTime()) - (xiFlight + lambdaFlight + negFlightPr);
      negDeltaTimeAsOmPi = (nTra.tofSignal() - nTra.tofEvTime()) - (omFlight + lambdaFlight + negFlightPi);
      negDeltaTimeAsOmPr = (nTra.tofSignal() - nTra.tofEvTime()) - (omFlight + lambdaFlight + negFlightPr);
    }
    if (bTra.hasTOF()) {
      bachDeltaTimeAsXiPi = (bTra.tofSignal() - bTra.tofEvTime()) - (xiFlight + bachFlightPi);
      bachDeltaTimeAsOmKa = (bTra.tofSignal() - bTra.tofEvTime()) - (omFlight + bachFlightKa);
    }

    casctofpids(
      posDeltaTimeAsXiPi, posDeltaTimeAsXiPr, negDeltaTimeAsXiPi, negDeltaTimeAsXiPr, bachDeltaTimeAsXiPi,
      posDeltaTimeAsOmPi, posDeltaTimeAsOmPr, negDeltaTimeAsOmPi, negDeltaTimeAsOmPr, bachDeltaTimeAsOmKa);

    float nSigmaXiLaPr = o2::aod::cascdata::kNoTOFValue;
    float nSigmaXiLaPi = o2::aod::cascdata::kNoTOFValue;
    float nSigmaXiPi = o2::aod::cascdata::kNoTOFValue;
    float nSigmaOmLaPr = o2::aod::cascdata::kNoTOFValue;
    float nSigmaOmLaPi = o2::aod::cascdata::kNoTOFValue;
    float nSigmaOmKa = o2::aod::cascdata::kNoTOFValue;

    // go for Nsigma values if requested
    if (doNSigmas && nSigmaCalibLoaded) {
      // Xi hypothesis ________________________
      if (cascade.sign() < 0) {         // XiMinus
        if (posDeltaTimeAsXiPr > -1e+5) // proton from Lambda from XiMinus has signal
          nSigmaXiLaPr = (posDeltaTimeAsXiPr - hMeanPosXiPr->Interpolate(cascade.p())) / hSigmaPosXiPr->Interpolate(cascade.p());
        if (negDeltaTimeAsXiPi > -1e+5) // pion from Lambda from XiMinus has signal
          nSigmaXiLaPi = (negDeltaTimeAsXiPi - hMeanNegXiPi->Interpolate(cascade.p())) / hSigmaNegXiPi->Interpolate(cascade.p());
        if (bachDeltaTimeAsXiPi > -1e+5) // pion from XiMinus has signal
          nSigmaXiPi = (bachDeltaTimeAsXiPi - hMeanBachXiPi->Interpolate(cascade.p())) / hSigmaBachXiPi->Interpolate(cascade.p());
        if (posDeltaTimeAsOmPr > -1e+5) // proton from Lambda from OmegaMinus has signal
          nSigmaOmLaPr = (posDeltaTimeAsOmPr - hMeanPosOmPr->Interpolate(cascade.p())) / hSigmaPosOmPr->Interpolate(cascade.p());
        if (negDeltaTimeAsOmPi > -1e+5) // pion from Lambda from OmegaMinus has signal
          nSigmaOmLaPi = (negDeltaTimeAsOmPi - hMeanNegOmPi->Interpolate(cascade.p())) / hSigmaNegOmPi->Interpolate(cascade.p());
        if (bachDeltaTimeAsOmKa > -1e+5) // kaon from OmegaMinus has signal
          nSigmaOmKa = (bachDeltaTimeAsOmKa - hMeanBachOmKa->Interpolate(cascade.p())) / hSigmaBachOmKa->Interpolate(cascade.p());
      } else {
        if (posDeltaTimeAsXiPi > -1e+5) // proton from Lambda from XiMinus has signal
          nSigmaXiLaPi = (posDeltaTimeAsXiPi - hMeanPosXiPi->Interpolate(cascade.p())) / hSigmaPosXiPi->Interpolate(cascade.p());
        if (negDeltaTimeAsXiPr > -1e+5) // pion from Lambda from XiMinus has signal
          nSigmaXiLaPr = (negDeltaTimeAsXiPr - hMeanNegXiPr->Interpolate(cascade.p())) / hSigmaNegXiPr->Interpolate(cascade.p());
        if (bachDeltaTimeAsXiPi > -1e+5) // pion from XiMinus has signal
          nSigmaXiPi = (bachDeltaTimeAsXiPi - hMeanBachXiPi->Interpolate(cascade.p())) / hSigmaBachXiPi->Interpolate(cascade.p());
        if (posDeltaTimeAsOmPi > -1e+5) // proton from Lambda from OmegaMinus has signal
          nSigmaOmLaPi = (posDeltaTimeAsOmPi - hMeanPosOmPi->Interpolate(cascade.p())) / hSigmaPosOmPi->Interpolate(cascade.p());
        if (negDeltaTimeAsOmPr > -1e+5) // pion from Lambda from OmegaMinus has signal
          nSigmaOmLaPr = (negDeltaTimeAsOmPr - hMeanNegOmPr->Interpolate(cascade.p())) / hSigmaNegOmPr->Interpolate(cascade.p());
        if (bachDeltaTimeAsOmKa > -1e+5) // kaon from OmegaMinus has signal
          nSigmaOmKa = (bachDeltaTimeAsOmKa - hMeanBachOmKa->Interpolate(cascade.p())) / hSigmaBachOmKa->Interpolate(cascade.p());
      }
      casctofnsigmas(nSigmaXiLaPi, nSigmaXiLaPr, nSigmaXiPi, nSigmaOmLaPi, nSigmaOmLaPr, nSigmaOmKa);
    }

    if (doQA) {
      // length factor due to eta (to offset e-loss)
      float positiveCosine = 1.0f / sqrt(1.0f + posTrack.getTgl() * posTrack.getTgl());
      float negativeCosine = 1.0f / sqrt(1.0f + negTrack.getTgl() * negTrack.getTgl());
      float bachelorCosine = 1.0f / sqrt(1.0f + bachTrack.getTgl() * bachTrack.getTgl());
      if (correctELossInclination.value == false) {
        negativeCosine = positiveCosine = bachelorCosine = 1.0f;
      }

      if (cascade.dcaV0daughters() < cascadeGroup.qaV0DCADau && cascade.dcacascdaughters() < cascadeGroup.qaCascDCADau && cascade.v0cosPA(collision.getX(), collision.getY(), collision.getZ()) > cascadeGroup.qaV0CosPA && cascade.casccosPA(collision.getX(), collision.getY(), collision.getZ()) > cascadeGroup.qaCascCosPA) {
        if (cascade.sign() < 0) {
          if (std::abs(cascade.mXi() - 1.32171) < cascadeGroup.qaMassWindow && fabs(pTra.tpcNSigmaPr()) < cascadeGroup.qaTPCNSigma && fabs(nTra.tpcNSigmaPi()) < cascadeGroup.qaTPCNSigma && fabs(bTra.tpcNSigmaPi()) < cascadeGroup.qaTPCNSigma && ((cascpdg == 0) || (cascpdg == 3312))) {
            histos.fill(HIST("h2dposDeltaTimeAsXiPr"), cascade.p(), cascade.eta(), posDeltaTimeAsXiPr);
            histos.fill(HIST("h2dnegDeltaTimeAsXiPi"), cascade.p(), cascade.eta(), negDeltaTimeAsXiPi);
            histos.fill(HIST("h2dbachDeltaTimeAsXiPi"), cascade.p(), cascade.eta(), bachDeltaTimeAsXiPi);
            if (calculationMethod.value == 2) {
              if (std::abs(posFlightPr_Method0 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon && std::abs(posFlightPr_Method1 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon) {
                histos.fill(HIST("hDeltaTimeMethodsVsP_posXiPr"), positiveP, cascade.positiveeta(), (posFlightPr_Method0 - posFlightPr_Method1) * positiveCosine);
                histos.fill(HIST("hRatioTimeMethodsVsP_posXiPr"), positiveP, cascade.positiveeta(), (posFlightPr_Method1 / posFlightPr_Method0) * positiveCosine);
              }
              if (std::abs(negFlightPi_Method0 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon && std::abs(negFlightPi_Method1 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon) {
                histos.fill(HIST("hDeltaTimeMethodsVsP_negXiPi"), negativeP, cascade.negativeeta(), (negFlightPi_Method0 - negFlightPi_Method1) * negativeCosine);
                histos.fill(HIST("hRatioTimeMethodsVsP_negXiPi"), negativeP, cascade.negativeeta(), (negFlightPi_Method1 / negFlightPi_Method0) * negativeCosine);
              }
              if (std::abs(bachFlightPi_Method0 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon && std::abs(bachFlightPi_Method1 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon) {
                histos.fill(HIST("hDeltaTimeMethodsVsP_bachXiPi"), bachelorP, cascade.bacheloreta(), (bachFlightPi_Method0 - bachFlightPi_Method1) * bachelorCosine);
                histos.fill(HIST("hRatioTimeMethodsVsP_bachXiPi"), bachelorP, cascade.bacheloreta(), (bachFlightPi_Method1 / bachFlightPi_Method0) * bachelorCosine);
              }
            }
            if (doQANSigma) {
              histos.fill(HIST("h2dNSigmaXiLaPi"), cascade.p(), nSigmaXiLaPi);
              histos.fill(HIST("h2dNSigmaXiLaPr"), cascade.p(), nSigmaXiLaPr);
              histos.fill(HIST("h2dNSigmaXiPi"), cascade.p(), nSigmaXiPi);
            }
          }
          if (std::abs(cascade.mOmega() - 1.67245) < cascadeGroup.qaMassWindow && fabs(pTra.tpcNSigmaPr()) < cascadeGroup.qaTPCNSigma && fabs(nTra.tpcNSigmaPi()) < cascadeGroup.qaTPCNSigma && fabs(bTra.tpcNSigmaKa()) < cascadeGroup.qaTPCNSigma && ((cascpdg == 0) || (cascpdg == 3334))) {
            histos.fill(HIST("h2dposDeltaTimeAsOmPr"), cascade.p(), cascade.eta(), posDeltaTimeAsOmPr);
            histos.fill(HIST("h2dnegDeltaTimeAsOmPi"), cascade.p(), cascade.eta(), negDeltaTimeAsOmPi);
            histos.fill(HIST("h2dbachDeltaTimeAsOmKa"), cascade.p(), cascade.eta(), bachDeltaTimeAsOmKa);
            if (calculationMethod.value == 2) {
              if (std::abs(posFlightPr_Method0 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon && std::abs(posFlightPr_Method1 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon) {
                histos.fill(HIST("hDeltaTimeMethodsVsP_posOmPr"), positiveP, cascade.positiveeta(), (posFlightPr_Method0 - posFlightPr_Method1) * positiveCosine);
                histos.fill(HIST("hRatioTimeMethodsVsP_posOmPr"), positiveP, cascade.positiveeta(), (posFlightPr_Method1 / posFlightPr_Method0) * positiveCosine);
              }
              if (std::abs(negFlightPi_Method0 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon && std::abs(negFlightPi_Method1 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon) {
                histos.fill(HIST("hDeltaTimeMethodsVsP_negOmPi"), negativeP, cascade.negativeeta(), (negFlightPi_Method0 - negFlightPi_Method1) * negativeCosine);
                histos.fill(HIST("hRatioTimeMethodsVsP_negOmPi"), negativeP, cascade.negativeeta(), (negFlightPi_Method1 / negFlightPi_Method0) * negativeCosine);
              }
              if (std::abs(bachFlightKa_Method0 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon && std::abs(bachFlightKa_Method1 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon) {
                histos.fill(HIST("hDeltaTimeMethodsVsP_bachOmKa"), bachelorP, cascade.bacheloreta(), (bachFlightKa_Method0 - bachFlightKa_Method1) * bachelorCosine);
                histos.fill(HIST("hRatioTimeMethodsVsP_bachOmKa"), bachelorP, cascade.bacheloreta(), (bachFlightKa_Method1 / bachFlightKa_Method0) * bachelorCosine);
              }
            }
            if (doQANSigma) {
              histos.fill(HIST("h2dNSigmaOmLaPi"), cascade.p(), nSigmaOmLaPi);
              histos.fill(HIST("h2dNSigmaOmLaPr"), cascade.p(), nSigmaOmLaPr);
              histos.fill(HIST("h2dNSigmaOmKa"), cascade.p(), nSigmaOmKa);
            }
          }
        } else {
          if (std::abs(cascade.mXi() - 1.32171) < cascadeGroup.qaMassWindow && fabs(pTra.tpcNSigmaPi()) < cascadeGroup.qaTPCNSigma && fabs(nTra.tpcNSigmaPr()) < cascadeGroup.qaTPCNSigma && fabs(bTra.tpcNSigmaPi()) < cascadeGroup.qaTPCNSigma && ((cascpdg == 0) || (cascpdg == -3312))) {
            histos.fill(HIST("h2dposDeltaTimeAsXiPi"), cascade.p(), cascade.eta(), posDeltaTimeAsXiPi);
            histos.fill(HIST("h2dnegDeltaTimeAsXiPr"), cascade.p(), cascade.eta(), negDeltaTimeAsXiPr);
            histos.fill(HIST("h2dbachDeltaTimeAsXiPi"), cascade.p(), cascade.eta(), bachDeltaTimeAsXiPi);
            if (calculationMethod.value == 2) {
              if (std::abs(posFlightPi_Method0 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon && std::abs(posFlightPi_Method1 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon) {
                histos.fill(HIST("hDeltaTimeMethodsVsP_posXiPi"), positiveP, cascade.positiveeta(), (posFlightPi_Method0 - posFlightPi_Method1) * positiveCosine);
                histos.fill(HIST("hRatioTimeMethodsVsP_posXiPi"), positiveP, cascade.positiveeta(), (posFlightPi_Method1 / posFlightPi_Method1) * positiveCosine);
              }
              if (std::abs(negFlightPr_Method0 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon && std::abs(negFlightPr_Method1 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon) {
                histos.fill(HIST("hDeltaTimeMethodsVsP_negXiPr"), negativeP, cascade.negativeeta(), (negFlightPr_Method0 - negFlightPr_Method1) * negativeCosine);
                histos.fill(HIST("hRatioTimeMethodsVsP_negXiPr"), negativeP, cascade.negativeeta(), (negFlightPr_Method1 / negFlightPr_Method0) * negativeCosine);
              }
              if (std::abs(bachFlightPi_Method0 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon && std::abs(bachFlightPi_Method1 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon) {
                histos.fill(HIST("hDeltaTimeMethodsVsP_bachXiPi"), bachelorP, cascade.bacheloreta(), (bachFlightPi_Method0 - bachFlightPi_Method1) * bachelorCosine);
                histos.fill(HIST("hRatioTimeMethodsVsP_bachXiPi"), bachelorP, cascade.bacheloreta(), (bachFlightPi_Method1 / bachFlightPi_Method0) * bachelorCosine);
              }
            }
            if (doQANSigma) {
              histos.fill(HIST("h2dNSigmaXiLaPi"), cascade.p(), nSigmaXiLaPi);
              histos.fill(HIST("h2dNSigmaXiLaPr"), cascade.p(), nSigmaXiLaPr);
              histos.fill(HIST("h2dNSigmaXiPi"), cascade.p(), nSigmaXiPi);
            }
          }
          if (std::abs(cascade.mOmega() - 1.67245) < cascadeGroup.qaMassWindow && fabs(pTra.tpcNSigmaPi()) < cascadeGroup.qaTPCNSigma && fabs(nTra.tpcNSigmaPr()) < cascadeGroup.qaTPCNSigma && fabs(bTra.tpcNSigmaKa()) < cascadeGroup.qaTPCNSigma && ((cascpdg == 0) || (cascpdg == -3334))) {
            histos.fill(HIST("h2dposDeltaTimeAsOmPi"), cascade.p(), cascade.eta(), posDeltaTimeAsOmPi);
            histos.fill(HIST("h2dnegDeltaTimeAsOmPr"), cascade.p(), cascade.eta(), negDeltaTimeAsOmPr);
            histos.fill(HIST("h2dbachDeltaTimeAsOmKa"), cascade.p(), cascade.eta(), bachDeltaTimeAsOmKa);
            if (calculationMethod.value == 2) {
              if (std::abs(posFlightPi_Method0 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon && std::abs(posFlightPi_Method1 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon) {
                histos.fill(HIST("hDeltaTimeMethodsVsP_posOmPi"), positiveP, cascade.positiveeta(), (posFlightPi_Method0 - posFlightPi_Method1) * positiveCosine);
                histos.fill(HIST("hRatioTimeMethodsVsP_posOmPi"), positiveP, cascade.positiveeta(), (posFlightPi_Method1 / posFlightPi_Method1) * positiveCosine);
              }
              if (std::abs(negFlightPr_Method0 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon && std::abs(negFlightPr_Method1 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon) {
                histos.fill(HIST("hDeltaTimeMethodsVsP_negOmPr"), negativeP, cascade.negativeeta(), (negFlightPr_Method0 - negFlightPr_Method1) * negativeCosine);
                histos.fill(HIST("hRatioTimeMethodsVsP_negOmPr"), negativeP, cascade.negativeeta(), (negFlightPr_Method1 / negFlightPr_Method0) * negativeCosine);
              }
              if (std::abs(bachFlightKa_Method0 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon && std::abs(bachFlightKa_Method1 - o2::aod::cascdata::kNoTOFValue) > o2::aod::cascdata::kEpsilon) {
                histos.fill(HIST("hDeltaTimeMethodsVsP_bachOmKa"), bachelorP, cascade.bacheloreta(), (bachFlightKa_Method0 - bachFlightKa_Method1) * bachelorCosine);
                histos.fill(HIST("hRatioTimeMethodsVsP_bachOmKa"), bachelorP, cascade.bacheloreta(), (bachFlightKa_Method1 / bachFlightKa_Method1) * bachelorCosine);
              }
            }
            if (doQANSigma) {
              histos.fill(HIST("h2dNSigmaOmLaPi"), cascade.p(), nSigmaOmLaPi);
              histos.fill(HIST("h2dNSigmaOmLaPr"), cascade.p(), nSigmaOmLaPr);
              histos.fill(HIST("h2dNSigmaOmKa"), cascade.p(), nSigmaOmKa);
            }
          }
        }
      }
    }
  }

  void processStandardData(aod::Collisions const& collisions, V0OriginalDatas const& V0s, CascOriginalDatas const& cascades, TracksWithAllExtras const&, aod::BCsWithTimestamps const& /*bcs*/)
  {
    // Fire up CCDB with first collision in record. If no collisions, bypass
    if (useCustomRunNumber || collisions.size() < 1) {
      initCCDB(manualRunNumber);
    } else {
      auto collision = collisions.begin();
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc.runNumber());
    }

    if (calculateV0s.value) {
      for (const auto& V0 : V0s) {
        // for storing whatever is the relevant quantity for the PV
        o2::dataformats::VertexBase primaryVertex;
        if (V0.has_collision()) {
          auto const& collision = V0.collision();
          primaryVertex.setPos({collision.posX(), collision.posY(), collision.posZ()});
          primaryVertex.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
        } else {
          primaryVertex.setPos({mVtx->getX(), mVtx->getY(), mVtx->getZ()});
        }

        auto pTra = V0.posTrack_as<TracksWithAllExtras>();
        auto nTra = V0.negTrack_as<TracksWithAllExtras>();
        processV0Candidate(primaryVertex, V0, pTra, nTra, 0);
      }
    }

    if (calculateCascades.value) {
      for (const auto& cascade : cascades) {
        // for storing whatever is the relevant quantity for the PV
        o2::dataformats::VertexBase primaryVertex;
        if (cascade.has_collision()) {
          auto const& collision = cascade.collision();
          primaryVertex.setPos({collision.posX(), collision.posY(), collision.posZ()});
          primaryVertex.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
        } else {
          primaryVertex.setPos({mVtx->getX(), mVtx->getY(), mVtx->getZ()});
        }

        auto pTra = cascade.posTrack_as<TracksWithAllExtras>();
        auto nTra = cascade.negTrack_as<TracksWithAllExtras>();
        auto bTra = cascade.bachelor_as<TracksWithAllExtras>();
        processCascadeCandidate(primaryVertex, cascade, pTra, nTra, bTra, 0);
      }
    }
  }

  void processDerivedData(soa::Join<aod::StraCollisions, aod::StraStamps> const& collisions, V0DerivedDatas const& V0s, CascDerivedDatas const& cascades, dauTracks const&)
  {
    // Fire up CCDB with first collision in record. If no collisions, bypass
    if (useCustomRunNumber || collisions.size() < 1) {
      initCCDB(manualRunNumber);
    } else {
      auto collision = collisions.begin();
      initCCDB(collision.runNumber());
    }

    if (calculateV0s.value) {
      for (const auto& V0 : V0s) {
        // for storing whatever is the relevant quantity for the PV
        o2::dataformats::VertexBase primaryVertex;
        if (V0.has_straCollision()) {
          auto const& collision = V0.straCollision_as<soa::Join<aod::StraCollisions, aod::StraStamps>>();
          primaryVertex.setPos({collision.posX(), collision.posY(), collision.posZ()});
          // cov: won't be used anyways, all fine
          primaryVertex.setCov(1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6);
        } else {
          primaryVertex.setPos({mVtx->getX(), mVtx->getY(), mVtx->getZ()});
        }

        auto pTra = V0.posTrackExtra_as<dauTracks>();
        auto nTra = V0.negTrackExtra_as<dauTracks>();
        processV0Candidate(primaryVertex, V0, pTra, nTra, 0);
      }
    }

    if (calculateCascades.value) {
      for (const auto& cascade : cascades) {
        // for storing whatever is the relevant quantity for the PV
        o2::dataformats::VertexBase primaryVertex;
        if (cascade.has_straCollision()) {
          auto const& collision = cascade.straCollision_as<soa::Join<aod::StraCollisions, aod::StraStamps>>();
          primaryVertex.setPos({collision.posX(), collision.posY(), collision.posZ()});
          primaryVertex.setCov(1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6);
        } else {
          primaryVertex.setPos({mVtx->getX(), mVtx->getY(), mVtx->getZ()});
        }

        auto pTra = cascade.posTrackExtra_as<dauTracks>();
        auto nTra = cascade.negTrackExtra_as<dauTracks>();
        auto bTra = cascade.bachTrackExtra_as<dauTracks>();
        processCascadeCandidate(primaryVertex, cascade, pTra, nTra, bTra, 0);
      }
    }
  }

  void processDerivedDataMCTest(soa::Join<aod::StraCollisions, aod::StraStamps> const& collisions, V0DerivedDatasMC const& V0s, CascDerivedDatasMC const& cascades, dauTracks const&, aod::V0MCCores const& v0mcs, aod::CascMCCores const& cascmcs)
  {
    // Fire up CCDB with first collision in record. If no collisions, bypass
    if (useCustomRunNumber || collisions.size() < 1) {
      initCCDB(manualRunNumber);
    } else {
      auto collision = collisions.begin();
      initCCDB(collision.runNumber());
    }

    if (calculateV0s.value) {
      for (const auto& V0 : V0s) {
        // for storing whatever is the relevant quantity for the PV
        o2::dataformats::VertexBase primaryVertex;
        if (V0.has_straCollision()) {
          auto const& collision = V0.straCollision_as<soa::Join<aod::StraCollisions, aod::StraStamps>>();
          primaryVertex.setPos({collision.posX(), collision.posY(), collision.posZ()});
          // cov: won't be used anyways, all fine
          primaryVertex.setCov(1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6);
        } else {
          primaryVertex.setPos({mVtx->getX(), mVtx->getY(), mVtx->getZ()});
        }

        // check association
        int v0pdg = 0;
        if (V0.v0MCCoreId() > -1) {
          auto v0mc = v0mcs.rawIteratorAt(V0.v0MCCoreId());
          v0pdg = v0mc.pdgCode();
          if (std::abs(v0pdg) != 3122 && v0pdg != 310) {
            continue; // only associated from this point on
          }
        }

        auto pTra = V0.posTrackExtra_as<dauTracks>();
        auto nTra = V0.negTrackExtra_as<dauTracks>();
        processV0Candidate(primaryVertex, V0, pTra, nTra, v0pdg);
      }
    }

    if (calculateCascades.value) {
      for (const auto& cascade : cascades) {
        // for storing whatever is the relevant quantity for the PV
        o2::dataformats::VertexBase primaryVertex;
        if (cascade.has_straCollision()) {
          auto const& collision = cascade.straCollision_as<soa::Join<aod::StraCollisions, aod::StraStamps>>();
          primaryVertex.setPos({collision.posX(), collision.posY(), collision.posZ()});
          primaryVertex.setCov(1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6);
        } else {
          primaryVertex.setPos({mVtx->getX(), mVtx->getY(), mVtx->getZ()});
        }

        // check association
        int cascpdg = 0;
        if (cascade.cascMCCoreId() > -1) {
          auto cascmc = cascmcs.rawIteratorAt(cascade.cascMCCoreId());
          cascpdg = cascmc.pdgCode();
          if (std::abs(cascpdg) != 3312 && std::abs(cascpdg) != 3334) {
            continue; // only associated from this point on
          }
        }

        auto pTra = cascade.posTrackExtra_as<dauTracks>();
        auto nTra = cascade.negTrackExtra_as<dauTracks>();
        auto bTra = cascade.bachTrackExtra_as<dauTracks>();
        processCascadeCandidate(primaryVertex, cascade, pTra, nTra, bTra, cascpdg);
      }
    }
  }

  PROCESS_SWITCH(strangenesstofpid, processStandardData, "Process standard data", false);
  PROCESS_SWITCH(strangenesstofpid, processDerivedData, "Process derived data", true);
  PROCESS_SWITCH(strangenesstofpid, processDerivedDataMCTest, "Process derived data / MC with assoc / not for analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenesstofpid>(cfgc)};
}
