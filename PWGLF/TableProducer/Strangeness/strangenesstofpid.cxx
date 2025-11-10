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
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
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
#include <string>
#include <utility>
#include <vector>

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
  Configurable<bool> reassociateTracks{"reassociateTracks", true, "if true, reassociate tracks to the collision the V0 or cascade belongs to. Relevant especially at high IR"};
  Configurable<bool> doBCshift{"doBCshift", true, "if true, perform time shift for collisions in different BCs when reassigning"};
  Configurable<bool> rejectUndefinedTof{"rejectUndefinedTof", true, "if true, reject tracks with TOF signal 0.000f for safety"};

  // auxiliary / debug tables as desired
  Configurable<int> calculateV0TOFPIDs{"calculateV0TOFPIDs", -1, "calculate V0TOFPIDs table (0: no, 1: yes, -1: auto)"};
  Configurable<int> calculateV0TOFBetas{"calculateV0TOFBetas", -1, "calculate V0TOFBetas table (0: no, 1: yes, -1: auto)"};
  Configurable<int> calculateV0TOFDebugs{"calculateV0TOFDebugs", -1, "calculate V0TOFDebugs table (0: no, 1: yes, -1: auto)"};
  Configurable<int> calculateCascTOFPIDs{"calculateCascTOFPIDs", -1, "calculate CascTOFPIDs table (0: no, 1: yes, -1: auto)"};

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

  struct : ConfigurableGroup {
    ConfigurableAxis axisPosition{"axisPosition", {400, -400.f, +400.f}, "position (cm)"};
    ConfigurableAxis axisEta{"axisEta", {20, -1.0f, +1.0f}, "#eta"};
    ConfigurableAxis axisDeltaTime{"axisDeltaTime", {2000, -1000.0f, +1000.0f}, "delta-time (ps)"};
    ConfigurableAxis axisDeltaTimeVsPrimaryCalculation{"axisDeltaTimeVsPrimaryCalculation", {500, -500.0f, +500.0f}, "delta-time (ps)"};
    ConfigurableAxis axisTime{"axisTime", {400, 10000.0f, +50000.0f}, "T (ps)"};
    ConfigurableAxis axisNSigma{"axisNSigma", {200, -10.0f, +10.0f}, "N(#sigma)"};
    ConfigurableAxis axisRatioMethods{"axisRatioMethods", {400, 0.9f, 1.9f}, "T_{method 1}/T_{method 0}"};
    ConfigurableAxis axisSnp{"axisSnp", {220, -1.1f, 1.1f}, "snp"};

    // master p axis
    ConfigurableAxis axisP{"axisP", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};

    // for zooming in at low values only (e-loss studies and effective correction)
    ConfigurableAxis axisSmallP{"axisSmallP", {250, 0.0f, 2.5f}, "p_{T} (GeV/c)"};

    // for BC shift QA plots
    // binning to actually match BC shifts but show in picoseconds
    const double bcShiftValuePS = o2::constants::lhc::LHCBunchSpacingNS * 1000.0f;
    ConfigurableAxis axisBCshift{"axisBCshift", {130, -120.5f * bcShiftValuePS, 9.5f * bcShiftValuePS}, "time shift (ps)"};

    // very broad time axis
    ConfigurableAxis axisTimeLong{"axisTimeLong", {3000, -1500000.0f, 1500000.0f}, "time (ps)"};
  } axes; // aggregate axes fo simplicity of navigation in HY

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

  // bookkeep propagation failures and successes
  enum typesOfPropagation { kPropagPosV0 = 0,
                            kPropagNegV0,
                            kPropagPosCasc,
                            kPropagNegCasc,
                            kPropagBachCasc,
                            kPropagTypes };

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
    if (calculateV0TOFPIDs.value < 0) {
      // check if TOF information is required, enable if so
      calculateV0TOFPIDs.value = isTableRequiredInWorkflow(initContext, "V0TOFPIDs");
      if (calculateV0TOFPIDs.value > 0) {
        LOGF(info, "Strangeness TOF PID: V0TOFPIDs calculations enabled automatically");
      }
    }
    if (calculateV0TOFBetas.value < 0) {
      // check if TOF information is required, enable if so
      calculateV0TOFBetas.value = isTableRequiredInWorkflow(initContext, "V0TOFBetas");
      if (calculateV0TOFBetas.value > 0) {
        LOGF(info, "Strangeness TOF PID: V0TOFBetas calculations enabled automatically");
      }
    }
    if (calculateV0TOFDebugs.value < 0) {
      // check if TOF information is required, enable if so
      calculateV0TOFDebugs.value = isTableRequiredInWorkflow(initContext, "V0TOFDebugs");
      if (calculateV0TOFDebugs.value > 0) {
        LOGF(info, "Strangeness TOF PID: V0TOFDebugs calculations enabled automatically");
      }
    }
    if (calculateCascTOFPIDs.value < 0) {
      // check if TOF information is required, enable if so
      calculateCascTOFPIDs.value = isTableRequiredInWorkflow(initContext, "CascTOFPIDs");
      if (calculateCascTOFPIDs.value > 0) {
        LOGF(info, "Strangeness TOF PID: CascTOFPIDs calculations enabled automatically");
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

    histos.add("hV0PositiveBCShift", "hV0PositiveBCShift", kTH1F, {axes.axisBCshift});
    histos.add("hV0NegativeBCShift", "hV0NegativeBCShift", kTH1F, {axes.axisBCshift});
    histos.add("hCascadePositiveBCShift", "hCascadePositiveBCShift", kTH1F, {axes.axisBCshift});
    histos.add("hCascadeNegativeBCShift", "hCascadeNegativeBCShift", kTH1F, {axes.axisBCshift});
    histos.add("hCascadeBachelorBCShift", "hCascadeBachelorBCShift", kTH1F, {axes.axisBCshift});

    histos.add("hTOFSignalPositive", "hTOFSignalPositive", kTH1F, {axes.axisTimeLong});
    histos.add("hTOFSignalNegative", "hTOFSignalNegative", kTH1F, {axes.axisTimeLong});

    histos.add("h2dTOFSignalPositive", "h2dTOFSignalPositive", kTH2F, {axes.axisTimeLong, axes.axisBCshift});
    histos.add("h2dTOFSignalNegative", "h2dTOFSignalNegative", kTH2F, {axes.axisTimeLong, axes.axisBCshift});

    histos.add("h2dTOFSignalCascadePositive", "h2dTOFSignalCascadePositive", kTH2F, {axes.axisTimeLong, axes.axisBCshift});
    histos.add("h2dTOFSignalCascadeNegative", "h2dTOFSignalCascadeNegative", kTH2F, {axes.axisTimeLong, axes.axisBCshift});
    histos.add("h2dTOFSignalCascadeBachelor", "h2dTOFSignalCascadeBachelor", kTH2F, {axes.axisTimeLong, axes.axisBCshift});

    histos.add("hCollisionTimes", "hCollisionTimes", kTH1F, {{2000, -1000.0f, 1000.0f}});

    // measured vs expected total time QA
    if (doQA) {
      // if in mode 1, bookkeep the failures of propagation
      if (calculationMethod.value == 1) {
        histos.add("hPropagationBookkeeping", "hPropagationBookkeeping", kTProfile, {{5, -0.5f, 4.5f}});
      }

      // standard deltaTime values
      if (calculateV0s.value > 0) {
        histos.add("h2dDeltaTimePositiveLambdaPi", "h2dDeltaTimePositiveLambdaPi", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dDeltaTimeNegativeLambdaPi", "h2dDeltaTimeNegativeLambdaPi", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dDeltaTimePositiveLambdaPr", "h2dDeltaTimePositiveLambdaPr", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dDeltaTimeNegativeLambdaPr", "h2dDeltaTimeNegativeLambdaPr", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dDeltaTimePositiveK0ShortPi", "h2dDeltaTimePositiveK0ShortPi", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dDeltaTimeNegativeK0ShortPi", "h2dDeltaTimeNegativeK0ShortPi", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});

        // delta time with respect to primary-like calculation
        histos.add("h2dDiffFromPrimCalcPositiveLambdaPi", "h2dDiffFromPrimCalcPositiveLambdaPi", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dDiffFromPrimCalcNegativeLambdaPi", "h2dDiffFromPrimCalcNegativeLambdaPi", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dDiffFromPrimCalcPositiveLambdaPr", "h2dDiffFromPrimCalcPositiveLambdaPr", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dDiffFromPrimCalcNegativeLambdaPr", "h2dDiffFromPrimCalcNegativeLambdaPr", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dDiffFromPrimCalcPositiveK0ShortPi", "h2dDiffFromPrimCalcPositiveK0ShortPi", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dDiffFromPrimCalcNegativeK0ShortPi", "h2dDiffFromPrimCalcNegativeK0ShortPi", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});

        // QA collision reassociation fraction (from track -> V0/cascade coll index)
        histos.add("h2dCorrectAssocPositiveLambdaPi", "h2dCorrectAssocPositiveLambdaPi", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dCorrectAssocNegativeLambdaPi", "h2dCorrectAssocNegativeLambdaPi", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dCorrectAssocPositiveLambdaPr", "h2dCorrectAssocPositiveLambdaPr", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dCorrectAssocNegativeLambdaPr", "h2dCorrectAssocNegativeLambdaPr", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dCorrectAssocPositiveK0ShortPi", "h2dCorrectAssocPositiveK0ShortPi", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dCorrectAssocNegativeK0ShortPi", "h2dCorrectAssocNegativeK0ShortPi", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
      }

      if (calculateCascades.value > 0) {
        histos.add("h2dposDeltaTimeAsXiPi", "h2dposDeltaTimeAsXiPi", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dposDeltaTimeAsXiPr", "h2dposDeltaTimeAsXiPr", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dnegDeltaTimeAsXiPi", "h2dnegDeltaTimeAsXiPi", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dnegDeltaTimeAsXiPr", "h2dnegDeltaTimeAsXiPr", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dbachDeltaTimeAsXiPi", "h2dbachDeltaTimeAsXiPi", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});

        histos.add("h2dposDeltaTimeAsOmPi", "h2dposDeltaTimeAsOmPi", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dposDeltaTimeAsOmPr", "h2dposDeltaTimeAsOmPr", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dnegDeltaTimeAsOmPi", "h2dnegDeltaTimeAsOmPi", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dnegDeltaTimeAsOmPr", "h2dnegDeltaTimeAsOmPr", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});
        histos.add("h2dbachDeltaTimeAsOmKa", "h2dbachDeltaTimeAsOmKa", {HistType::kTH3F, {axes.axisP, axes.axisEta, axes.axisDeltaTime}});

        // delta time with respect to primary-like calculation
        histos.add("h2dposDiffFromPrimCalcAsXiPi", "h2dposDiffFromPrimCalcAsXiPi", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dposDiffFromPrimCalcAsXiPr", "h2dposDiffFromPrimCalcAsXiPr", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dnegDiffFromPrimCalcAsXiPi", "h2dnegDiffFromPrimCalcAsXiPi", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dnegDiffFromPrimCalcAsXiPr", "h2dnegDiffFromPrimCalcAsXiPr", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dbachDiffFromPrimCalcAsXiPi", "h2dbachDiffFromPrimCalcAsXiPi", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});

        histos.add("h2dposDiffFromPrimCalcAsOmPi", "h2dposDiffFromPrimCalcAsOmPi", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dposDiffFromPrimCalcAsOmPr", "h2dposDiffFromPrimCalcAsOmPr", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dnegDiffFromPrimCalcAsOmPi", "h2dnegDiffFromPrimCalcAsOmPi", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dnegDiffFromPrimCalcAsOmPr", "h2dnegDiffFromPrimCalcAsOmPr", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});
        histos.add("h2dbachDiffFromPrimCalcAsOmKa", "h2dbachDiffFromPrimCalcAsOmKa", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTimeVsPrimaryCalculation}});

        // QA collision reassociation fraction (from track -> V0/cascade coll index)
        histos.add("h2dposCorrectAssocAsXiPi", "h2dposCorrectAssocAsXiPi", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dposCorrectAssocAsXiPr", "h2dposCorrectAssocAsXiPr", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dnegCorrectAssocAsXiPi", "h2dnegCorrectAssocAsXiPi", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dnegCorrectAssocAsXiPr", "h2dnegCorrectAssocAsXiPr", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dbachCorrectAssocAsXiPi", "h2dbachCorrectAssocAsXiPi", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});

        histos.add("h2dposCorrectAssocAsOmPi", "h2dposCorrectAssocAsOmPi", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dposCorrectAssocAsOmPr", "h2dposCorrectAssocAsOmPr", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dnegCorrectAssocAsOmPi", "h2dnegCorrectAssocAsOmPi", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dnegCorrectAssocAsOmPr", "h2dnegCorrectAssocAsOmPr", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
        histos.add("h2dbachCorrectAssocAsOmKa", "h2dbachCorrectAssocAsOmKa", {HistType::kTH2F, {axes.axisP, {2, -0.5f, 1.5f}}});
      }

      histos.add("h2dPositiveTOFProperties", "h2dPositiveTOFProperties", {HistType::kTH2F, {axes.axisP, {4, -0.5, 3.5f}}});
      histos.add("h2dNegativeTOFProperties", "h2dNegativeTOFProperties", {HistType::kTH2F, {axes.axisP, {4, -0.5, 3.5f}}});

      if (doQANSigma) {
        if (calculateV0s.value > 0) {
          histos.add("h2dNSigmaPositiveLambdaPi", "h2dNSigmaPositiveLambdaPi", {HistType::kTH2F, {axes.axisP, axes.axisNSigma}});
          histos.add("h2dNSigmaNegativeLambdaPi", "h2dNSigmaNegativeLambdaPi", {HistType::kTH2F, {axes.axisP, axes.axisNSigma}});
          histos.add("h2dNSigmaPositiveLambdaPr", "h2dNSigmaPositiveLambdaPr", {HistType::kTH2F, {axes.axisP, axes.axisNSigma}});
          histos.add("h2dNSigmaNegativeLambdaPr", "h2dNSigmaNegativeLambdaPr", {HistType::kTH2F, {axes.axisP, axes.axisNSigma}});
          histos.add("h2dNSigmaPositiveK0ShortPi", "h2dNSigmaPositiveK0ShortPi", {HistType::kTH2F, {axes.axisP, axes.axisNSigma}});
          histos.add("h2dNSigmaNegativeK0ShortPi", "h2dNSigmaNegativeK0ShortPi", {HistType::kTH2F, {axes.axisP, axes.axisNSigma}});
        }

        if (calculateCascades.value > 0) {
          histos.add("h2dNSigmaXiLaPi", "h2dNSigmaXiLaPi", {HistType::kTH2F, {axes.axisP, axes.axisNSigma}});
          histos.add("h2dNSigmaXiLaPr", "h2dNSigmaXiLaPr", {HistType::kTH2F, {axes.axisP, axes.axisNSigma}});
          histos.add("h2dNSigmaXiPi", "h2dNSigmaXiPi", {HistType::kTH2F, {axes.axisP, axes.axisNSigma}});
          histos.add("h2dNSigmaOmLaPi", "h2dNSigmaOmLaPi", {HistType::kTH2F, {axes.axisP, axes.axisNSigma}});
          histos.add("h2dNSigmaOmLaPr", "h2dNSigmaOmLaPr", {HistType::kTH2F, {axes.axisP, axes.axisNSigma}});
          histos.add("h2dNSigmaOmKa", "h2dNSigmaOmKa", {HistType::kTH2F, {axes.axisP, axes.axisNSigma}});
        }
      }

      // delta lambda decay time
      histos.add("h2dLambdaDeltaDecayTime", "h2dLambdaDeltaDecayTime", {HistType::kTH2F, {axes.axisP, axes.axisDeltaTime}});
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
    float timeK0Short = o2::aod::v0data::kNoTOFValue;
    float timeLambda = o2::aod::v0data::kNoTOFValue;
    float timePositivePr = o2::aod::v0data::kNoTOFValue;
    float timePositivePi = o2::aod::v0data::kNoTOFValue;
    float timeNegativePr = o2::aod::v0data::kNoTOFValue;
    float timeNegativePi = o2::aod::v0data::kNoTOFValue;

    float timeAsPrimaryPositivePr = o2::aod::v0data::kNoTOFValue;
    float timeAsPrimaryPositivePi = o2::aod::v0data::kNoTOFValue;
    float timeAsPrimaryNegativePr = o2::aod::v0data::kNoTOFValue;
    float timeAsPrimaryNegativePi = o2::aod::v0data::kNoTOFValue;

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

    // extra auxiliary variables
    float deltaDecayTimeLambda = o2::aod::v0data::kNoTOFValue;
    float deltaDecayTimeAntiLambda = o2::aod::v0data::kNoTOFValue;
    float deltaDecayTimeK0Short = o2::aod::v0data::kNoTOFValue;

    float betaLambda = o2::aod::v0data::kNoTOFValue;
    float betaAntiLambda = o2::aod::v0data::kNoTOFValue;
    float betaK0Short = o2::aod::v0data::kNoTOFValue;
  };

  // structs to hold information
  struct cascTofInfo { // holds processed information regarding TOF for Cascades
    float posFlightPi = o2::aod::cascdata::kNoTOFValue;
    float posFlightPr = o2::aod::cascdata::kNoTOFValue;
    float negFlightPi = o2::aod::cascdata::kNoTOFValue;
    float negFlightPr = o2::aod::cascdata::kNoTOFValue;
    float bachFlightPi = o2::aod::cascdata::kNoTOFValue;
    float bachFlightKa = o2::aod::cascdata::kNoTOFValue;

    float posFlightAsPrimaryPi = o2::aod::cascdata::kNoTOFValue;
    float posFlightAsPrimaryPr = o2::aod::cascdata::kNoTOFValue;
    float negFlightAsPrimaryPi = o2::aod::cascdata::kNoTOFValue;
    float negFlightAsPrimaryPr = o2::aod::cascdata::kNoTOFValue;
    float bachFlightAsPrimaryPi = o2::aod::cascdata::kNoTOFValue;
    float bachFlightAsPrimaryKa = o2::aod::cascdata::kNoTOFValue;

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
    v0tof.timeK0Short = lengthV0 / velocityK0Short; // in picoseconds
    v0tof.timeLambda = lengthV0 / velocityLambda;   // in picoseconds

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

    bool pValidTOF = rejectUndefinedTof.value ? static_cast<bool>(std::fabs(pTof.tofSignal) > o2::aod::v0data::kEpsilon) : true;
    bool nValidTOF = rejectUndefinedTof.value ? static_cast<bool>(std::fabs(nTof.tofSignal) > o2::aod::v0data::kEpsilon) : true;

    //_____________________________________________________________________________________________
    // Actual calculation
    float velocityPositivePr, velocityPositivePi, lengthPositive;
    velocityPositivePr = velocityPositivePi = lengthPositive = o2::aod::v0data::kNoTOFValue;

    if (pTof.hasTOF && pTof.tofEvTime > -1e+5 && pValidTOF) {
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
          if (doQA) {
            histos.fill(HIST("hPropagationBookkeeping"), kPropagPosV0, static_cast<float>(successPropag));
          }
          if (successPropag) {
            lengthPositive = pTof.length - ltIntegral.getL();
            v0tof.timePositivePr = o2::framework::pid::tof::MassToExpTime(pTof.tofExpMom, lengthPositive, o2::constants::physics::MassProton * o2::constants::physics::MassProton);
            v0tof.timePositivePi = o2::framework::pid::tof::MassToExpTime(pTof.tofExpMom, lengthPositive, o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);

            // as primary
            v0tof.timeAsPrimaryPositivePr = o2::framework::pid::tof::MassToExpTime(pTof.tofExpMom, pTof.length, o2::constants::physics::MassProton * o2::constants::physics::MassProton);
            v0tof.timeAsPrimaryPositivePi = o2::framework::pid::tof::MassToExpTime(pTof.tofExpMom, pTof.length, o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);
          }
        }
      }
      if (lengthPositive > 0.0f) {
        v0tof.deltaTimePositiveLambdaPr = (pTof.tofSignal - pTof.tofEvTime) - (v0tof.timeLambda + v0tof.timePositivePr);
        v0tof.deltaTimePositiveLambdaPi = (pTof.tofSignal - pTof.tofEvTime) - (v0tof.timeLambda + v0tof.timePositivePi);
        v0tof.deltaTimePositiveK0ShortPi = (pTof.tofSignal - pTof.tofEvTime) - (v0tof.timeK0Short + v0tof.timePositivePi);

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
              histos.fill(HIST("h2dCorrectAssocPositiveLambdaPr"), v0.p(), static_cast<float>(collisionId == pTof.collisionId));
              histos.fill(HIST("h2dDiffFromPrimCalcPositiveLambdaPr"), v0.p(), (pTof.tofSignal - pTof.tofEvTime) - v0tof.timeAsPrimaryPositivePr);
              if (doQANSigma && std::fabs(v0tof.nSigmaPositiveLambdaPr - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaPositiveLambdaPr"), v0.pt(), v0tof.nSigmaPositiveLambdaPr);
              }
            }
            if (antiLambdaCandidate) {
              histos.fill(HIST("h2dDeltaTimePositiveLambdaPi"), v0.p(), v0.eta(), v0tof.deltaTimePositiveLambdaPi);
              histos.fill(HIST("h2dCorrectAssocPositiveLambdaPi"), v0.p(), static_cast<float>(collisionId == pTof.collisionId));
              histos.fill(HIST("h2dDiffFromPrimCalcPositiveLambdaPi"), v0.p(), (pTof.tofSignal - pTof.tofEvTime) - v0tof.timeAsPrimaryPositivePi);
              if (doQANSigma && std::fabs(v0tof.nSigmaPositiveLambdaPi - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaPositiveLambdaPi"), v0.pt(), v0tof.nSigmaPositiveLambdaPi);
              }
            }
            if (k0ShortCandidate) {
              histos.fill(HIST("h2dDeltaTimePositiveK0ShortPi"), v0.p(), v0.eta(), v0tof.deltaTimePositiveK0ShortPi);
              histos.fill(HIST("h2dCorrectAssocPositiveK0ShortPi"), v0.p(), static_cast<float>(collisionId == pTof.collisionId));
              histos.fill(HIST("h2dDiffFromPrimCalcPositiveK0ShortPi"), v0.p(), (pTof.tofSignal - pTof.tofEvTime) - v0tof.timeAsPrimaryPositivePi);
              if (doQANSigma && std::fabs(v0tof.nSigmaPositiveK0ShortPi - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaPositiveK0ShortPi"), v0.pt(), v0tof.nSigmaPositiveK0ShortPi);
              }
            }
          }
        }
      }
    }
    float velocityNegativePr, velocityNegativePi, lengthNegative;
    velocityNegativePr = velocityNegativePi = lengthNegative = o2::aod::v0data::kNoTOFValue;
    if (nTof.hasTOF && nTof.tofEvTime > -1e+5 && nValidTOF) {
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
          if (doQA) {
            histos.fill(HIST("hPropagationBookkeeping"), kPropagNegV0, static_cast<float>(successPropag));
          }
          if (successPropag) {
            lengthNegative = nTof.length - ltIntegral.getL();
            v0tof.timeNegativePr = o2::framework::pid::tof::MassToExpTime(nTof.tofExpMom, lengthNegative, o2::constants::physics::MassProton * o2::constants::physics::MassProton);
            v0tof.timeNegativePi = o2::framework::pid::tof::MassToExpTime(nTof.tofExpMom, lengthNegative, o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);

            // as primary
            v0tof.timeAsPrimaryNegativePr = o2::framework::pid::tof::MassToExpTime(nTof.tofExpMom, nTof.length, o2::constants::physics::MassProton * o2::constants::physics::MassProton);
            v0tof.timeAsPrimaryNegativePi = o2::framework::pid::tof::MassToExpTime(nTof.tofExpMom, nTof.length, o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);
          }
        }
      }
      if (lengthNegative > 0.0f) {
        v0tof.deltaTimeNegativeLambdaPr = (nTof.tofSignal - nTof.tofEvTime) - (v0tof.timeLambda + v0tof.timeNegativePr);
        v0tof.deltaTimeNegativeLambdaPi = (nTof.tofSignal - nTof.tofEvTime) - (v0tof.timeLambda + v0tof.timeNegativePi);
        v0tof.deltaTimeNegativeK0ShortPi = (nTof.tofSignal - nTof.tofEvTime) - (v0tof.timeK0Short + v0tof.timeNegativePi);

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
              histos.fill(HIST("h2dCorrectAssocNegativeLambdaPi"), v0.p(), static_cast<float>(collisionId == nTof.collisionId));
              histos.fill(HIST("h2dDiffFromPrimCalcNegativeLambdaPi"), v0.p(), (nTof.tofSignal - nTof.tofEvTime) - v0tof.timeAsPrimaryNegativePi);
              if (doQANSigma && std::fabs(v0tof.nSigmaNegativeLambdaPi - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaNegativeLambdaPi"), v0.pt(), v0tof.nSigmaNegativeLambdaPi);
              }
            }
            if (antiLambdaCandidate) {
              histos.fill(HIST("h2dDeltaTimeNegativeLambdaPr"), v0.p(), v0.eta(), v0tof.deltaTimeNegativeLambdaPr);
              histos.fill(HIST("h2dCorrectAssocNegativeLambdaPr"), v0.p(), static_cast<float>(collisionId == nTof.collisionId));
              histos.fill(HIST("h2dDiffFromPrimCalcNegativeLambdaPr"), v0.p(), (nTof.tofSignal - nTof.tofEvTime) - v0tof.timeAsPrimaryNegativePr);
              if (doQANSigma && std::fabs(v0tof.nSigmaNegativeLambdaPr - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaNegativeLambdaPr"), v0.pt(), v0tof.nSigmaNegativeLambdaPr);
              }
            }
            if (k0ShortCandidate) {
              histos.fill(HIST("h2dDeltaTimeNegativeK0ShortPi"), v0.p(), v0.eta(), v0tof.deltaTimeNegativeK0ShortPi);
              histos.fill(HIST("h2dCorrectAssocNegativeK0ShortPi"), v0.p(), static_cast<float>(collisionId == nTof.collisionId));
              histos.fill(HIST("h2dDiffFromPrimCalcNegativeK0ShortPi"), v0.p(), (nTof.tofSignal - nTof.tofEvTime) - v0tof.timeAsPrimaryNegativePi);
              if (doQANSigma && std::fabs(v0tof.nSigmaNegativeK0ShortPi - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaNegativeK0ShortPi"), v0.pt(), v0tof.nSigmaNegativeK0ShortPi);
              }
            }
          }
        }
      }

      // temporarily commented out
      // bool compatibleK0Short = true;
      // int incompatibilityReason = 0;
      // if (std::abs(v0tof.nSigmaPositiveK0ShortPi - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon && std::abs(v0tof.nSigmaPositiveK0ShortPi) > 4) {
      //   compatibleK0Short = false; // reject only if info present and incompatible
      //   incompatibilityReason += 1;
      // }
      // if (std::abs(v0tof.nSigmaNegativeK0ShortPi - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon && std::abs(v0tof.nSigmaNegativeK0ShortPi) > 4) {
      //   compatibleK0Short = false; // reject only if info present and incompatible
      //   incompatibilityReason += 2;
      // }

      // if(!compatibleK0Short && passesQAcuts && k0ShortCandidate){
      //   histos.fill(HIST("hIncompatibilityReason"), incompatibilityReason);
      //   // LOGF(info, "Incompatible K0, sigmas = (%.2f %.2f), lengths = (%.2f %.2f) tofSignals = (%.2f %.2f) evtimes = (%.2f %.2f)", v0tof.nSigmaPositiveK0ShortPi, v0tof.nSigmaNegativeK0ShortPi, lengthPositive, lengthNegative, pTof.tofSignal, nTof.tofSignal, pTof.tofEvTime, nTof.tofEvTime);
      // }

      // calculation of delta-decay-time (no reliance on event time)
      if (nTof.hasTOF && pTof.hasTOF > 0) { // does not depend on event time
        v0tof.deltaDecayTimeLambda = (pTof.tofSignal - v0tof.timePositivePr) - (nTof.tofSignal - v0tof.timeNegativePi);
        v0tof.deltaDecayTimeAntiLambda = (pTof.tofSignal - v0tof.timePositivePi) - (nTof.tofSignal - v0tof.timeNegativePr);
        v0tof.deltaDecayTimeK0Short = (pTof.tofSignal - v0tof.timePositivePi) - (nTof.tofSignal - v0tof.timeNegativePi);

        float evTimeMean = 0.5f * (pTof.tofEvTime + nTof.tofEvTime);
        float decayTimeLambda = 0.5f * ((pTof.tofSignal - v0tof.timePositivePr) + (nTof.tofSignal - v0tof.timeNegativePi)) - evTimeMean;
        float decayTimeAntiLambda = 0.5f * ((pTof.tofSignal - v0tof.timePositivePi) + (nTof.tofSignal - v0tof.timeNegativePr)) - evTimeMean;
        float decayTimeK0Short = 0.5f * ((pTof.tofSignal - v0tof.timePositivePi) + (nTof.tofSignal - v0tof.timeNegativePi)) - evTimeMean;

        constexpr float lightSpeed = 0.0299792458; // in cm/ps
        v0tof.betaLambda = (lengthV0 / decayTimeLambda) / lightSpeed;
        v0tof.betaAntiLambda = (lengthV0 / decayTimeAntiLambda) / lightSpeed;
        v0tof.betaK0Short = (lengthV0 / decayTimeK0Short) / lightSpeed;
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

    bool pValidTOF = rejectUndefinedTof.value ? static_cast<bool>(std::fabs(pTof.tofSignal) > o2::aod::v0data::kEpsilon) : true;
    bool nValidTOF = rejectUndefinedTof.value ? static_cast<bool>(std::fabs(nTof.tofSignal) > o2::aod::v0data::kEpsilon) : true;
    bool bValidTOF = rejectUndefinedTof.value ? static_cast<bool>(std::fabs(bTof.tofSignal) > o2::aod::v0data::kEpsilon) : true;

    //_____________________________________________________________________________________________
    // Actual calculation
    if (pTof.hasTOF && pTof.tofEvTime > -1e+5 && pValidTOF) {
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
          if (doQA) {
            histos.fill(HIST("hPropagationBookkeeping"), kPropagPosCasc, static_cast<float>(successPropag));
          }
          if (successPropag) {
            lengthPositive = pTof.length - ltIntegral.getL();
            casctof.posFlightPr = o2::framework::pid::tof::MassToExpTime(pTof.tofExpMom, pTof.length - ltIntegral.getL(), o2::constants::physics::MassProton * o2::constants::physics::MassProton);
            casctof.posFlightPi = o2::framework::pid::tof::MassToExpTime(pTof.tofExpMom, pTof.length - ltIntegral.getL(), o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);

            // as primary
            casctof.posFlightAsPrimaryPr = o2::framework::pid::tof::MassToExpTime(pTof.tofExpMom, pTof.length, o2::constants::physics::MassProton * o2::constants::physics::MassProton);
            casctof.posFlightAsPrimaryPi = o2::framework::pid::tof::MassToExpTime(pTof.tofExpMom, pTof.length, o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);
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
              histos.fill(HIST("h2dposCorrectAssocAsXiPr"), cascade.p(), static_cast<float>(collisionId == pTof.collisionId));
              histos.fill(HIST("h2dposDiffFromPrimCalcAsXiPr"), cascade.p(), (pTof.tofSignal - pTof.tofEvTime) - casctof.posFlightAsPrimaryPr);
              if (doQANSigma && std::fabs(casctof.nSigmaXiLaPr - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaXiLaPr"), cascade.pt(), casctof.nSigmaXiLaPr);
              }
            }
            if (xiPlusCandidate) {
              histos.fill(HIST("h2dposDeltaTimeAsXiPi"), cascade.p(), cascade.eta(), casctof.posDeltaTimeAsXiPi);
              histos.fill(HIST("h2dposCorrectAssocAsXiPi"), cascade.p(), static_cast<float>(collisionId == pTof.collisionId));
              histos.fill(HIST("h2dposDiffFromPrimCalcAsXiPi"), cascade.p(), (pTof.tofSignal - pTof.tofEvTime) - casctof.posFlightAsPrimaryPi);
              if (doQANSigma && std::fabs(casctof.nSigmaXiLaPi - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaXiLaPi"), cascade.pt(), casctof.nSigmaXiLaPi);
              }
            }
            if (omegaMinusCandidate) {
              histos.fill(HIST("h2dposDeltaTimeAsOmPr"), cascade.p(), cascade.eta(), casctof.posDeltaTimeAsOmPr);
              histos.fill(HIST("h2dposCorrectAssocAsOmPr"), cascade.p(), static_cast<float>(collisionId == pTof.collisionId));
              histos.fill(HIST("h2dposDiffFromPrimCalcAsOmPr"), cascade.p(), (pTof.tofSignal - pTof.tofEvTime) - casctof.posFlightAsPrimaryPr);
              if (doQANSigma && std::fabs(casctof.nSigmaOmLaPr - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaOmLaPr"), cascade.pt(), casctof.nSigmaOmLaPr);
              }
            }
            if (omegaPlusCandidate) {
              histos.fill(HIST("h2dposDeltaTimeAsOmPi"), cascade.p(), cascade.eta(), casctof.posDeltaTimeAsOmPi);
              histos.fill(HIST("h2dposCorrectAssocAsOmPi"), cascade.p(), static_cast<float>(collisionId == pTof.collisionId));
              histos.fill(HIST("h2dposDiffFromPrimCalcAsOmPi"), cascade.p(), (pTof.tofSignal - pTof.tofEvTime) - casctof.posFlightAsPrimaryPi);
              if (doQANSigma && std::fabs(casctof.nSigmaOmLaPi - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaOmLaPi"), cascade.pt(), casctof.nSigmaOmLaPi);
              }
            }
          }
        }
      }
    } // end positive

    if (nTof.hasTOF && nTof.tofEvTime > -1e+5 && nValidTOF) {
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
          if (doQA) {
            histos.fill(HIST("hPropagationBookkeeping"), kPropagNegCasc, static_cast<float>(successPropag));
          }
          if (successPropag) {
            lengthNegative = nTof.length - ltIntegral.getL();
            casctof.negFlightPr = o2::framework::pid::tof::MassToExpTime(nTof.tofExpMom, nTof.length - ltIntegral.getL(), o2::constants::physics::MassProton * o2::constants::physics::MassProton);
            casctof.negFlightPi = o2::framework::pid::tof::MassToExpTime(nTof.tofExpMom, nTof.length - ltIntegral.getL(), o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);

            // as primary
            casctof.negFlightAsPrimaryPr = o2::framework::pid::tof::MassToExpTime(nTof.tofExpMom, nTof.length, o2::constants::physics::MassProton * o2::constants::physics::MassProton);
            casctof.negFlightAsPrimaryPi = o2::framework::pid::tof::MassToExpTime(nTof.tofExpMom, nTof.length, o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);
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
            casctof.nSigmaXiLaPi = (casctof.negDeltaTimeAsXiPi - hMeanNegXiPi->Interpolate(cascade.p())) / hSigmaNegXiPi->Interpolate(cascade.p());
            casctof.nSigmaOmLaPi = (casctof.negDeltaTimeAsOmPi - hMeanNegOmPi->Interpolate(cascade.p())) / hSigmaNegOmPi->Interpolate(cascade.p());
          } else {
            casctof.nSigmaXiLaPr = (casctof.negDeltaTimeAsXiPr - hMeanNegXiPr->Interpolate(cascade.p())) / hSigmaNegXiPr->Interpolate(cascade.p());
            casctof.nSigmaOmLaPr = (casctof.negDeltaTimeAsOmPr - hMeanNegOmPr->Interpolate(cascade.p())) / hSigmaNegOmPr->Interpolate(cascade.p());
          }
        }

        // do QA histograms (calibration / QC)
        if (doQA) {
          if (passesQAcuts) {
            if (xiMinusCandidate) {
              histos.fill(HIST("h2dnegDeltaTimeAsXiPi"), cascade.p(), cascade.eta(), casctof.negDeltaTimeAsXiPi);
              histos.fill(HIST("h2dnegCorrectAssocAsXiPi"), cascade.p(), static_cast<float>(collisionId == nTof.collisionId));
              histos.fill(HIST("h2dnegDiffFromPrimCalcAsXiPi"), cascade.p(), (nTof.tofSignal - nTof.tofEvTime) - casctof.negFlightAsPrimaryPi);
              if (doQANSigma && std::fabs(casctof.nSigmaXiLaPi - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaXiLaPi"), cascade.pt(), casctof.nSigmaXiLaPi);
              }
            }
            if (xiPlusCandidate) {
              histos.fill(HIST("h2dnegDeltaTimeAsXiPr"), cascade.p(), cascade.eta(), casctof.negDeltaTimeAsXiPr);
              histos.fill(HIST("h2dnegCorrectAssocAsXiPr"), cascade.p(), static_cast<float>(collisionId == nTof.collisionId));
              histos.fill(HIST("h2dnegDiffFromPrimCalcAsXiPr"), cascade.p(), (nTof.tofSignal - nTof.tofEvTime) - casctof.negFlightAsPrimaryPr);
              if (doQANSigma && std::fabs(casctof.nSigmaXiLaPr - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaXiLaPr"), cascade.pt(), casctof.nSigmaXiLaPr);
              }
            }
            if (omegaMinusCandidate) {
              histos.fill(HIST("h2dnegDeltaTimeAsOmPi"), cascade.p(), cascade.eta(), casctof.negDeltaTimeAsOmPi);
              histos.fill(HIST("h2dnegCorrectAssocAsOmPi"), cascade.p(), static_cast<float>(collisionId == nTof.collisionId));
              histos.fill(HIST("h2dnegDiffFromPrimCalcAsOmPi"), cascade.p(), (nTof.tofSignal - nTof.tofEvTime) - casctof.negFlightAsPrimaryPi);
              if (doQANSigma && std::fabs(casctof.nSigmaOmLaPi - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaOmLaPi"), cascade.pt(), casctof.nSigmaOmLaPi);
              }
            }
            if (omegaPlusCandidate) {
              histos.fill(HIST("h2dnegDeltaTimeAsOmPr"), cascade.p(), cascade.eta(), casctof.negDeltaTimeAsOmPr);
              histos.fill(HIST("h2dnegCorrectAssocAsOmPr"), cascade.p(), static_cast<float>(collisionId == nTof.collisionId));
              histos.fill(HIST("h2dnegDiffFromPrimCalcAsOmPr"), cascade.p(), (nTof.tofSignal - nTof.tofEvTime) - casctof.negFlightAsPrimaryPr);
              if (doQANSigma && std::fabs(casctof.nSigmaOmLaPr - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaOmLaPr"), cascade.pt(), casctof.nSigmaOmLaPr);
              }
            }
          }
        }
      }
    } // end negative

    if (bTof.hasTOF && bTof.tofEvTime > -1e+5 && bValidTOF) {
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
          if (doQA) {
            histos.fill(HIST("hPropagationBookkeeping"), kPropagBachCasc, static_cast<float>(successPropag));
          }
          if (successPropag) {
            lengthBachelor = bTof.length - ltIntegral.getL();
            casctof.bachFlightPi = o2::framework::pid::tof::MassToExpTime(bTof.tofExpMom, bTof.length - ltIntegral.getL(), o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);
            casctof.bachFlightKa = o2::framework::pid::tof::MassToExpTime(bTof.tofExpMom, bTof.length - ltIntegral.getL(), o2::constants::physics::MassKaonCharged * o2::constants::physics::MassKaonCharged);

            // as primary
            casctof.bachFlightAsPrimaryPi = o2::framework::pid::tof::MassToExpTime(bTof.tofExpMom, bTof.length, o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged);
            casctof.bachFlightAsPrimaryKa = o2::framework::pid::tof::MassToExpTime(bTof.tofExpMom, bTof.length, o2::constants::physics::MassKaonCharged * o2::constants::physics::MassKaonCharged);
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
              histos.fill(HIST("h2dbachCorrectAssocAsXiPi"), cascade.p(), static_cast<float>(collisionId == bTof.collisionId));
              histos.fill(HIST("h2dbachDiffFromPrimCalcAsXiPi"), cascade.p(), (bTof.tofSignal - bTof.tofEvTime) - casctof.bachFlightAsPrimaryPi);
              if (doQANSigma && std::fabs(casctof.nSigmaXiPi - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaXiPi"), cascade.pt(), casctof.nSigmaXiPi);
              }
            }
            if (xiPlusCandidate) {
              histos.fill(HIST("h2dbachDeltaTimeAsXiPi"), cascade.p(), cascade.eta(), casctof.bachDeltaTimeAsXiPi);
              histos.fill(HIST("h2dbachCorrectAssocAsXiPi"), cascade.p(), static_cast<float>(collisionId == bTof.collisionId));
              histos.fill(HIST("h2dbachDiffFromPrimCalcAsXiPi"), cascade.p(), (bTof.tofSignal - bTof.tofEvTime) - casctof.bachFlightAsPrimaryPi);
              if (doQANSigma && std::fabs(casctof.nSigmaXiPi - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaXiPi"), cascade.pt(), casctof.nSigmaXiPi);
              }
            }
            if (omegaMinusCandidate) {
              histos.fill(HIST("h2dbachDeltaTimeAsOmKa"), cascade.p(), cascade.eta(), casctof.bachDeltaTimeAsOmKa);
              histos.fill(HIST("h2dbachCorrectAssocAsOmKa"), cascade.p(), static_cast<float>(collisionId == bTof.collisionId));
              histos.fill(HIST("h2dbachDiffFromPrimCalcAsOmKa"), cascade.p(), (bTof.tofSignal - bTof.tofEvTime) - casctof.bachFlightAsPrimaryKa);
              if (doQANSigma && std::fabs(casctof.nSigmaOmKa - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaOmKa"), cascade.pt(), casctof.nSigmaOmKa);
              }
            }
            if (omegaPlusCandidate) {
              histos.fill(HIST("h2dbachDeltaTimeAsOmKa"), cascade.p(), cascade.eta(), casctof.bachDeltaTimeAsOmKa);
              histos.fill(HIST("h2dbachCorrectAssocAsOmKa"), cascade.p(), static_cast<float>(collisionId == bTof.collisionId));
              histos.fill(HIST("h2dbachDiffFromPrimCalcAsOmKa"), cascade.p(), (bTof.tofSignal - bTof.tofEvTime) - casctof.bachFlightAsPrimaryKa);
              if (doQANSigma && std::fabs(casctof.nSigmaOmKa - o2::aod::v0data::kNoTOFValue) > o2::aod::v0data::kEpsilon) {
                histos.fill(HIST("h2dNSigmaOmKa"), cascade.pt(), casctof.nSigmaOmKa);
              }
            }
          }
        }
      }
    } // end bachelor

    // don't forget to give feedback
    return casctof;
  }

  void processStandardData(aod::BCs const& bcs, aod::Collisions const& collisions, V0OriginalDatas const& V0s, CascOriginalDatas const& cascades, TracksWithAllExtras const& tracks, aod::BCsWithTimestamps const& /*bcs*/)
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
      if (track.hasTOF() && track.has_collision()) {
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
      histos.fill(HIST("hCollisionTimes"), collisionEventTime[collision.globalIndex()]);
    }

    if (calculateV0s.value) {
      for (const auto& V0 : V0s) {
        trackTofInfo pTof, nTof; // information storage

        auto pTra = V0.posTrack_as<TracksWithAllExtras>();
        auto nTra = V0.negTrack_as<TracksWithAllExtras>();
        double deltaTimePos = 0.0f;
        double deltaTimeNeg = 0.0f;

        auto collisionV0 = collisions.rawIteratorAt(V0.collisionId());
        auto bcV0 = bcs.rawIteratorAt(collisionV0.bcId());

        if (pTra.collisionId() >= 0) {
          auto collisionPos = collisions.rawIteratorAt(pTra.collisionId());
          auto bcPos = bcs.rawIteratorAt(collisionPos.bcId());
          const int64_t deltaBcPos = bcPos.globalBC() - bcV0.globalBC();
          deltaTimePos = o2::constants::lhc::LHCBunchSpacingNS * deltaBcPos * 1000.0f;
          histos.fill(HIST("hV0PositiveBCShift"), deltaTimePos);
        }

        if (nTra.collisionId() >= 0) {
          auto collisionNeg = collisions.rawIteratorAt(nTra.collisionId());
          auto bcNeg = bcs.rawIteratorAt(collisionNeg.bcId());
          const int64_t deltaBcNeg = bcNeg.globalBC() - bcV0.globalBC();
          deltaTimeNeg = o2::constants::lhc::LHCBunchSpacingNS * deltaBcNeg * 1000.0f;
          histos.fill(HIST("hV0NegativeBCShift"), deltaTimeNeg);
        }

        pTof.collisionId = pTra.collisionId();
        pTof.hasITS = pTra.hasITS();
        pTof.hasTPC = pTra.hasTPC();
        pTof.hasTOF = pTra.hasTOF();
        pTof.tofExpMom = pTra.tofExpMom();
        pTof.tofEvTime = reassociateTracks ? collisionEventTime[V0.collisionId()] : pTra.tofEvTime();
        pTof.tofSignal = pTra.tofSignal() + (doBCshift ? deltaTimePos : 0.0f);
        pTof.length = pTra.length();
        pTof.tpcNSigmaPi = pTra.tpcNSigmaPi();
        pTof.tpcNSigmaPr = pTra.tpcNSigmaPr();

        nTof.collisionId = nTra.collisionId();
        nTof.hasITS = nTra.hasITS();
        nTof.hasTPC = nTra.hasTPC();
        nTof.hasTOF = nTra.hasTOF();
        nTof.tofExpMom = nTra.tofExpMom();
        nTof.tofEvTime = reassociateTracks ? collisionEventTime[V0.collisionId()] : nTra.tofEvTime();
        nTof.tofSignal = nTra.tofSignal() + (doBCshift ? deltaTimeNeg : 0.0f);
        nTof.length = nTra.length();
        nTof.tpcNSigmaPi = nTra.tpcNSigmaPi();
        nTof.tpcNSigmaPr = nTra.tpcNSigmaPr();

        if (pTof.hasTOF) {
          histos.fill(HIST("hTOFSignalPositive"), pTof.tofSignal);
          histos.fill(HIST("h2dTOFSignalPositive"), pTof.tofSignal, deltaTimePos);
        }

        if (nTof.hasTOF) {
          histos.fill(HIST("hTOFSignalNegative"), nTof.tofSignal);
          histos.fill(HIST("h2dTOFSignalNegative"), nTof.tofSignal, deltaTimeNeg);
        }

        v0TofInfo v0tof = calculateTofInfoV0(collisions, V0.collisionId(), V0, pTof, nTof);

        if (doNSigmas) {
          v0tofnsigmas(
            v0tof.nSigmaPositiveLambdaPr, v0tof.nSigmaNegativeLambdaPi,
            v0tof.nSigmaNegativeLambdaPr, v0tof.nSigmaPositiveLambdaPi,
            v0tof.nSigmaPositiveK0ShortPi, v0tof.nSigmaNegativeK0ShortPi);
        }
        if (calculateV0TOFPIDs.value) {
          v0tofpid(v0tof.deltaTimePositiveLambdaPi, v0tof.deltaTimePositiveLambdaPr,
                   v0tof.deltaTimeNegativeLambdaPi, v0tof.deltaTimeNegativeLambdaPr,
                   v0tof.deltaTimePositiveK0ShortPi, v0tof.deltaTimeNegativeK0ShortPi,
                   v0tof.deltaDecayTimeLambda, v0tof.deltaDecayTimeAntiLambda, v0tof.deltaDecayTimeK0Short);
        }
        if (calculateV0TOFBetas.value) {
          v0tofbeta(v0tof.betaLambda, v0tof.betaAntiLambda, v0tof.betaK0Short);
        }
        if (calculateV0TOFDebugs.value) {
          v0tofdebugs(v0tof.timeLambda, v0tof.timeK0Short,
                      v0tof.timePositivePr, v0tof.timePositivePi,
                      v0tof.timeNegativePr, v0tof.timeNegativePi);
        }
      }
    }

    if (calculateCascades.value) {
      for (const auto& cascade : cascades) {
        trackTofInfo pTof, nTof, bTof; // information storage

        auto pTra = cascade.posTrack_as<TracksWithAllExtras>();
        auto nTra = cascade.negTrack_as<TracksWithAllExtras>();
        auto bTra = cascade.bachelor_as<TracksWithAllExtras>();

        double deltaTimePos = 0.0f;
        double deltaTimeNeg = 0.0f;
        double deltaTimeBach = 0.0f;

        auto collisionCascade = collisions.rawIteratorAt(cascade.collisionId());
        auto bcV0 = bcs.rawIteratorAt(collisionCascade.bcId());

        if (pTra.collisionId() >= 0) {
          auto collisionPos = collisions.rawIteratorAt(pTra.collisionId());
          auto bcPos = bcs.rawIteratorAt(collisionPos.bcId());
          const int64_t deltaBcPos = bcPos.globalBC() - bcV0.globalBC();
          deltaTimePos = o2::constants::lhc::LHCBunchSpacingNS * deltaBcPos * 1000.0f;
          histos.fill(HIST("hCascadePositiveBCShift"), deltaTimePos);
        }

        if (nTra.collisionId() >= 0) {
          auto collisionNeg = collisions.rawIteratorAt(nTra.collisionId());
          auto bcNeg = bcs.rawIteratorAt(collisionNeg.bcId());
          const int64_t deltaBcNeg = bcNeg.globalBC() - bcV0.globalBC();
          deltaTimeNeg = o2::constants::lhc::LHCBunchSpacingNS * deltaBcNeg * 1000.0f;
          histos.fill(HIST("hCascadeNegativeBCShift"), deltaTimeNeg);
        }

        if (bTra.collisionId() >= 0) {
          auto collisionBach = collisions.rawIteratorAt(bTra.collisionId());
          auto bcBach = bcs.rawIteratorAt(collisionBach.bcId());
          const int64_t deltaBcBach = bcBach.globalBC() - bcV0.globalBC();
          deltaTimeBach = o2::constants::lhc::LHCBunchSpacingNS * deltaBcBach * 1000.0f;
          histos.fill(HIST("hCascadeBachelorBCShift"), deltaTimeBach);
        }

        pTof.collisionId = pTra.collisionId();
        pTof.hasITS = pTra.hasITS();
        pTof.hasTPC = pTra.hasTPC();
        pTof.hasTOF = pTra.hasTOF();
        pTof.tofExpMom = pTra.tofExpMom();
        pTof.tofEvTime = reassociateTracks ? collisionEventTime[cascade.collisionId()] : pTra.tofEvTime();
        pTof.tofSignal = pTra.tofSignal() + (doBCshift ? deltaTimePos : 0.0f);
        pTof.length = pTra.length();
        pTof.tpcNSigmaPi = pTra.tpcNSigmaPi();
        pTof.tpcNSigmaPr = pTra.tpcNSigmaPr();

        nTof.collisionId = nTra.collisionId();
        nTof.hasITS = nTra.hasITS();
        nTof.hasTPC = nTra.hasTPC();
        nTof.hasTOF = nTra.hasTOF();
        nTof.tofExpMom = nTra.tofExpMom();
        nTof.tofEvTime = reassociateTracks ? collisionEventTime[cascade.collisionId()] : nTra.tofEvTime();
        nTof.tofSignal = nTra.tofSignal() + (doBCshift ? deltaTimeNeg : 0.0f);
        nTof.length = nTra.length();
        nTof.tpcNSigmaPi = nTra.tpcNSigmaPi();
        nTof.tpcNSigmaPr = nTra.tpcNSigmaPr();

        bTof.collisionId = bTra.collisionId();
        bTof.hasITS = bTra.hasITS();
        bTof.hasTPC = bTra.hasTPC();
        bTof.hasTOF = bTra.hasTOF();
        bTof.tofExpMom = bTra.tofExpMom();
        bTof.tofEvTime = reassociateTracks ? collisionEventTime[cascade.collisionId()] : bTra.tofEvTime();
        bTof.tofSignal = bTra.tofSignal() + (doBCshift ? deltaTimeBach : 0.0f);
        bTof.length = bTra.length();
        bTof.tpcNSigmaPi = bTra.tpcNSigmaPi();
        bTof.tpcNSigmaKa = bTra.tpcNSigmaKa();

        if (pTof.hasTOF) {
          histos.fill(HIST("h2dTOFSignalCascadePositive"), pTof.tofSignal, deltaTimePos);
        }
        if (nTof.hasTOF) {
          histos.fill(HIST("h2dTOFSignalCascadeNegative"), nTof.tofSignal, deltaTimeNeg);
        }
        if (bTof.hasTOF) {
          histos.fill(HIST("h2dTOFSignalCascadeBachelor"), bTof.tofSignal, deltaTimeBach);
        }

        cascTofInfo casctof = calculateTofInfoCascade(collisions, cascade.collisionId(), cascade, pTof, nTof, bTof);

        if (doNSigmas) {
          casctofnsigmas(
            casctof.nSigmaXiLaPi, casctof.nSigmaXiLaPr, casctof.nSigmaXiPi,
            casctof.nSigmaOmLaPi, casctof.nSigmaOmLaPr, casctof.nSigmaOmKa);
        }
        if (calculateCascTOFPIDs.value) {
          casctofpids(
            casctof.posDeltaTimeAsXiPi, casctof.posDeltaTimeAsXiPr,
            casctof.negDeltaTimeAsXiPi, casctof.negDeltaTimeAsXiPr, casctof.bachDeltaTimeAsXiPi,
            casctof.posDeltaTimeAsOmPi, casctof.posDeltaTimeAsOmPr,
            casctof.negDeltaTimeAsOmPi, casctof.negDeltaTimeAsOmPr, casctof.bachDeltaTimeAsOmKa);
        }
      }
    }
  }

  void processDerivedData(soa::Join<aod::StraCollisions, aod::StraStamps, aod::StraEvTimes> const& collisions, V0DerivedDatas const& V0s, CascDerivedDatas const& cascades, dauTracks const& dauTrackTable, aod::DauTrackTOFPIDs const& dauTrackTOFPIDs)
  {
    bool isNewTOFFormat = true; // can only happen for new format

    for (const auto& collision : collisions) {
      histos.fill(HIST("hCollisionTimes"), collision.eventTime());
    }

    // auto-determine if using old format
    if (dauTrackTOFPIDs.size() != 0) {
      auto firstTOFPID = dauTrackTOFPIDs.rawIteratorAt(0);
      isNewTOFFormat = firstTOFPID.straCollisionId() < 0 ? false : true;
    }

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

        double deltaTimeBcPos = 1e+6;
        double deltaTimeBcNeg = 1e+6;

        pTof.hasITS = pTra.hasITS();
        pTof.hasTPC = pTra.hasTPC();
        pTof.hasTOF = pTra.hasTOF();
        pTof.tpcNSigmaPi = pTra.tpcNSigmaPi();
        pTof.tpcNSigmaPr = pTra.tpcNSigmaPr();
        if (tofIndices[V0.posTrackExtraId()] >= 0 && collision.eventTime() > -1e+5) {
          auto pTofExt = dauTrackTOFPIDs.rawIteratorAt(tofIndices[V0.posTrackExtraId()]);

          if (pTofExt.straCollisionId() >= 0) {
            // extract BC for BC time shift
            auto collisionTrack = collisions.rawIteratorAt(pTofExt.straCollisionId());
            const int64_t deltaBc = collisionTrack.globalBC() - collision.globalBC();
            const double deltaTimeBc = o2::constants::lhc::LHCBunchSpacingNS * deltaBc * 1000.0f;
            histos.fill(HIST("hV0PositiveBCShift"), deltaTimeBc);
            deltaTimeBcPos = deltaTimeBc;

            // assign variables
            pTof.collisionId = pTofExt.straCollisionId();
            pTof.tofExpMom = pTofExt.tofExpMom();
            pTof.tofEvTime = reassociateTracks.value ? collision.eventTime() : pTofExt.tofEvTime();
            pTof.tofSignal = pTofExt.tofSignal() + (doBCshift.value ? deltaTimeBc : 0.0f);
            pTof.length = pTofExt.length();
          }
        }

        nTof.hasITS = nTra.hasITS();
        nTof.hasTPC = nTra.hasTPC();
        nTof.hasTOF = nTra.hasTOF();
        nTof.tpcNSigmaPi = nTra.tpcNSigmaPi();
        nTof.tpcNSigmaPr = nTra.tpcNSigmaPr();
        if (tofIndices[V0.negTrackExtraId()] >= 0 && collision.eventTime() > -1e+5) {
          auto nTofExt = dauTrackTOFPIDs.rawIteratorAt(tofIndices[V0.negTrackExtraId()]);

          if (nTofExt.straCollisionId() >= 0) {
            // extract BC for BC time shift
            auto collisionTrack = collisions.rawIteratorAt(nTofExt.straCollisionId());
            const int64_t deltaBc = collisionTrack.globalBC() - collision.globalBC();
            const double deltaTimeBc = o2::constants::lhc::LHCBunchSpacingNS * deltaBc * 1000.0f;
            histos.fill(HIST("hV0NegativeBCShift"), deltaTimeBc);
            deltaTimeBcNeg = deltaTimeBc;

            // assign variables
            nTof.collisionId = nTofExt.straCollisionId();
            nTof.tofExpMom = nTofExt.tofExpMom();
            nTof.tofEvTime = reassociateTracks.value ? collision.eventTime() : nTofExt.tofEvTime();
            nTof.tofSignal = nTofExt.tofSignal() + (doBCshift.value ? deltaTimeBc : 0.0f);
            nTof.length = nTofExt.length();
          }
        }
        if (pTof.hasTOF) {
          histos.fill(HIST("hTOFSignalPositive"), pTof.tofSignal);
          histos.fill(HIST("h2dTOFSignalPositive"), pTof.tofSignal, deltaTimeBcPos);
        }

        if (nTof.hasTOF) {
          histos.fill(HIST("hTOFSignalNegative"), nTof.tofSignal);
          histos.fill(HIST("h2dTOFSignalNegative"), nTof.tofSignal, deltaTimeBcNeg);
        }

        v0TofInfo v0tof = calculateTofInfoV0(collisions, V0.straCollisionId(), V0, pTof, nTof);

        if (doNSigmas) {
          v0tofnsigmas(
            v0tof.nSigmaPositiveLambdaPr, v0tof.nSigmaNegativeLambdaPi,
            v0tof.nSigmaNegativeLambdaPr, v0tof.nSigmaPositiveLambdaPi,
            v0tof.nSigmaPositiveK0ShortPi, v0tof.nSigmaNegativeK0ShortPi);
        }
        if (calculateV0TOFPIDs.value) {
          v0tofpid(v0tof.deltaTimePositiveLambdaPi, v0tof.deltaTimePositiveLambdaPr,
                   v0tof.deltaTimeNegativeLambdaPi, v0tof.deltaTimeNegativeLambdaPr,
                   v0tof.deltaTimePositiveK0ShortPi, v0tof.deltaTimeNegativeK0ShortPi,
                   v0tof.deltaDecayTimeLambda, v0tof.deltaDecayTimeAntiLambda, v0tof.deltaDecayTimeK0Short);
        }
        if (calculateV0TOFBetas.value) {
          v0tofbeta(v0tof.betaLambda, v0tof.betaAntiLambda, v0tof.betaK0Short);
        }
        if (calculateV0TOFDebugs.value) {
          v0tofdebugs(v0tof.timeLambda, v0tof.timeK0Short,
                      v0tof.timePositivePr, v0tof.timePositivePi,
                      v0tof.timeNegativePr, v0tof.timeNegativePi);
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

          if (pTofExt.straCollisionId() >= 0) {
            // extract BC for BC time shift
            auto collisionTrack = collisions.rawIteratorAt(pTofExt.straCollisionId());
            const int64_t deltaBc = collisionTrack.globalBC() - collision.globalBC();
            const double deltaTimeBc = o2::constants::lhc::LHCBunchSpacingNS * deltaBc * 1000.0f;
            histos.fill(HIST("hCascadePositiveBCShift"), deltaTimeBc);
            histos.fill(HIST("h2dTOFSignalCascadePositive"), pTof.tofSignal, deltaTimeBc);

            pTof.collisionId = pTofExt.straCollisionId();
            pTof.tofExpMom = pTofExt.tofExpMom();
            pTof.tofEvTime = reassociateTracks.value ? collision.eventTime() : pTofExt.tofEvTime();
            pTof.tofSignal = pTofExt.tofSignal() + (doBCshift.value ? deltaTimeBc : 0.0f);
            pTof.length = pTofExt.length();
          }
        }

        nTof.hasITS = nTra.hasITS();
        nTof.hasTPC = nTra.hasTPC();
        nTof.hasTOF = nTra.hasTOF();
        nTof.tpcNSigmaPi = nTra.tpcNSigmaPi();
        nTof.tpcNSigmaPr = nTra.tpcNSigmaPr();
        if (tofIndices[cascade.negTrackExtraId()] >= 0 && collision.eventTime() > -1e+5) {
          auto nTofExt = dauTrackTOFPIDs.rawIteratorAt(tofIndices[cascade.negTrackExtraId()]);

          if (nTofExt.straCollisionId() >= 0) {
            // extract BC for BC time shift
            auto collisionTrack = collisions.rawIteratorAt(nTofExt.straCollisionId());
            const int64_t deltaBc = collisionTrack.globalBC() - collision.globalBC();
            const double deltaTimeBc = o2::constants::lhc::LHCBunchSpacingNS * deltaBc * 1000.0f;
            histos.fill(HIST("hCascadeNegativeBCShift"), deltaTimeBc);
            histos.fill(HIST("h2dTOFSignalCascadeNegative"), nTof.tofSignal, deltaTimeBc);

            nTof.collisionId = nTofExt.straCollisionId();
            nTof.tofExpMom = nTofExt.tofExpMom();
            nTof.tofEvTime = reassociateTracks.value ? collision.eventTime() : nTofExt.tofEvTime();
            nTof.tofSignal = nTofExt.tofSignal() + (doBCshift.value ? deltaTimeBc : 0.0f);
            nTof.length = nTofExt.length();
          }
        }

        bTof.hasITS = bTra.hasITS();
        bTof.hasTPC = bTra.hasTPC();
        bTof.hasTOF = bTra.hasTOF();
        bTof.tpcNSigmaPi = bTra.tpcNSigmaPi();
        bTof.tpcNSigmaKa = bTra.tpcNSigmaKa();
        if (tofIndices[cascade.bachTrackExtraId()] >= 0 && collision.eventTime() > -1e+5) {
          auto bTofExt = dauTrackTOFPIDs.rawIteratorAt(tofIndices[cascade.bachTrackExtraId()]);

          if (bTofExt.straCollisionId() >= 0) {
            // extract BC for BC time shift
            auto collisionTrack = collisions.rawIteratorAt(bTofExt.straCollisionId());
            const int64_t deltaBc = collisionTrack.globalBC() - collision.globalBC();
            const double deltaTimeBc = o2::constants::lhc::LHCBunchSpacingNS * deltaBc * 1000.0f;
            histos.fill(HIST("hCascadeBachelorBCShift"), deltaTimeBc);
            histos.fill(HIST("h2dTOFSignalCascadeBachelor"), bTof.tofSignal, deltaTimeBc);

            bTof.collisionId = bTofExt.straCollisionId();
            bTof.tofExpMom = bTofExt.tofExpMom();
            bTof.tofEvTime = reassociateTracks.value ? collision.eventTime() : bTofExt.tofEvTime();
            bTof.tofSignal = bTofExt.tofSignal() + (doBCshift.value ? deltaTimeBc : 0.0f);
            bTof.length = bTofExt.length();
          }
        }

        cascTofInfo casctof = calculateTofInfoCascade(collisions, cascade.straCollisionId(), cascade, pTof, nTof, bTof);

        if (doNSigmas) {
          casctofnsigmas(
            casctof.nSigmaXiLaPi, casctof.nSigmaXiLaPr, casctof.nSigmaXiPi,
            casctof.nSigmaOmLaPi, casctof.nSigmaOmLaPr, casctof.nSigmaOmKa);
        }
        if (calculateCascTOFPIDs.value) {
          casctofpids(
            casctof.posDeltaTimeAsXiPi, casctof.posDeltaTimeAsXiPr,
            casctof.negDeltaTimeAsXiPi, casctof.negDeltaTimeAsXiPr, casctof.bachDeltaTimeAsXiPi,
            casctof.posDeltaTimeAsOmPi, casctof.posDeltaTimeAsOmPr,
            casctof.negDeltaTimeAsOmPi, casctof.negDeltaTimeAsOmPr, casctof.bachDeltaTimeAsOmKa);
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
