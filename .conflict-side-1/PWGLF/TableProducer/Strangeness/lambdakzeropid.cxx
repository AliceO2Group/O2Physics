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

#include <cmath>
#include <array>
#include <cstdlib>
#include <map>
#include <iterator>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// For original data loops
using V0OriginalDatas = soa::Join<aod::V0Indices, aod::V0Cores>;
using TracksWithAllExtras = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullHe, aod::TOFEvTime, aod::TOFSignal>;

// For derived data analysis
using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs, aod::DauTrackTOFPIDs>;
using V0DerivedDatas = soa::Join<aod::V0Cores, aod::V0Extras, aod::V0CollRefs>;

struct lambdakzeropid {
  // TOF pid for strangeness (recalculated with topology)
  Produces<aod::V0TOFPIDs> v0tofpid;        // table with Nsigmas
  Produces<aod::V0TOFBetas> v0tofbeta;      // table with betas
  Produces<aod::V0TOFDebugs> v0tofdebugs;   // table with extra debug information
  Produces<aod::V0TOFNSigmas> v0tofnsigmas; // table with nsigmas

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // mean vertex position to be used if no collision associated
  o2::dataformats::MeanVertexObject* mVtx = nullptr;

  // For manual sliceBy
  Preslice<V0OriginalDatas> perCollisionOriginal = o2::aod::v0data::collisionId;
  ;
  Preslice<V0DerivedDatas> perCollisionDerived = o2::aod::v0data::straCollisionId;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Operation and minimisation criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<float> tofPosition{"tofPosition", 377.934f, "TOF effective (inscribed) radius"};
  Configurable<bool> doQA{"doQA", true, "create QA histos"};
  Configurable<bool> doQANSigma{"doQANSigma", true, "create QA of Nsigma histos"};
  Configurable<float> qaDCADau{"qaDCADau", 0.5, "DCA daughters (cm) for QA plots"};
  Configurable<float> qaCosPA{"qaCosPA", 0.999, "CosPA for QA plots"};
  Configurable<float> qaMassWindow{"qaMassWindow", 0.005, "Mass window around expected (in GeV/c2) for QA plots"};
  Configurable<float> qaTPCNSigma{"qaTPCNSigma", 5, "TPC N-sigma to apply for qa plots"};
  Configurable<bool> doNSigmas{"doNSigmas", false, "calculate TOF N-sigma"};

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> nSigmaPath{"nSigmaPath", "Users/d/ddobrigk/stratof", "Path of information for n-sigma calculation"};
  Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

  // manual
  Configurable<int> useCustomRunNumber{"useCustomRunNumber", false, "Use custom timestamp"};
  Configurable<int> manualRunNumber{"manualRunNumber", 544122, "manual run number if no collisions saved"};

  ConfigurableAxis axisEta{"axisEta", {20, -1.0f, +1.0f}, "#eta"};
  ConfigurableAxis axisDeltaTime{"axisDeltaTime", {2000, -1000.0f, +1000.0f}, "delta-time (ps)"};
  ConfigurableAxis axisTime{"axisTime", {200, 0.0f, +20000.0f}, "T (ps)"};
  ConfigurableAxis axisNSigma{"axisNSigma", {200, -10.0f, +10.0f}, "N(#sigma)"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "p_{T} (GeV/c)"};

  // for n-sigma calibration
  bool nSigmaCalibLoaded;
  TList* nSigmaCalibObjects;
  TH1 *hMeanPosLaPi, *hSigmaPosLaPi;
  TH1 *hMeanPosLaPr, *hSigmaPosLaPr;
  TH1 *hMeanNegLaPi, *hSigmaNegLaPi;
  TH1 *hMeanNegLaPr, *hSigmaNegLaPr;
  TH1 *hMeanPosK0Pi, *hSigmaPosK0Pi;
  TH1 *hMeanNegK0Pi, *hSigmaNegK0Pi;

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

  void init(InitContext&)
  {
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

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // per event
    histos.add("hCandidateCounter", "hCandidateCounter", kTH1F, {{500, -0.5f, 499.5f}});

    // measured vs expected total time QA
    if (doQA) {
      histos.add("h2dProtonMeasuredVsExpected", "h2dProtonMeasuredVsExpected", {HistType::kTH2F, {axisTime, axisTime}});
      histos.add("h2dPionMeasuredVsExpected", "h2dPionMeasuredVsExpected", {HistType::kTH2F, {axisTime, axisTime}});

      // standard deltaTime values
      histos.add("h2dDeltaTimePositiveLambdaPi", "h2dDeltaTimePositiveLambdaPi", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dDeltaTimeNegativeLambdaPi", "h2dDeltaTimeNegativeLambdaPi", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dDeltaTimePositiveLambdaPr", "h2dDeltaTimePositiveLambdaPr", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dDeltaTimeNegativeLambdaPr", "h2dDeltaTimeNegativeLambdaPr", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dDeltaTimePositiveK0ShortPi", "h2dDeltaTimePositiveK0ShortPi", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});
      histos.add("h2dDeltaTimeNegativeK0ShortPi", "h2dDeltaTimeNegativeK0ShortPi", {HistType::kTH3F, {axisPt, axisEta, axisDeltaTime}});

      histos.add("h2dPositiveTOFProperties", "h2dPositiveTOFProperties", {HistType::kTH2F, {axisPt, {4, -0.5, 3.5f}}});
      histos.add("h2dNegativeTOFProperties", "h2dNegativeTOFProperties", {HistType::kTH2F, {axisPt, {4, -0.5, 3.5f}}});

      if (doQANSigma) {
        // standard NSigma values
        histos.add("h2dNSigmaPositiveLambdaPi", "h2dNSigmaPositiveLambdaPi", {HistType::kTH2F, {axisPt, axisNSigma}});
        histos.add("h2dNSigmaNegativeLambdaPi", "h2dNSigmaNegativeLambdaPi", {HistType::kTH2F, {axisPt, axisNSigma}});
        histos.add("h2dNSigmaPositiveLambdaPr", "h2dNSigmaPositiveLambdaPr", {HistType::kTH2F, {axisPt, axisNSigma}});
        histos.add("h2dNSigmaNegativeLambdaPr", "h2dNSigmaNegativeLambdaPr", {HistType::kTH2F, {axisPt, axisNSigma}});
        histos.add("h2dNSigmaPositiveK0ShortPi", "h2dNSigmaPositiveK0ShortPi", {HistType::kTH2F, {axisPt, axisNSigma}});
        histos.add("h2dNSigmaNegativeK0ShortPi", "h2dNSigmaNegativeK0ShortPi", {HistType::kTH2F, {axisPt, axisNSigma}});
      }

      // delta lambda decay time
      histos.add("h2dLambdaDeltaDecayTime", "h2dLambdaDeltaDecayTime", {HistType::kTH2F, {axisPt, axisDeltaTime}});
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
      mVtx = ccdb->getForRun<o2::dataformats::MeanVertexObject>(mVtxPath, runNumber);
      mRunNumber = runNumber;
      return;
    }

    o2::parameters::GRPObject* grpo = ccdb->getForRun<o2::parameters::GRPObject>(grpPath, runNumber);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for run " << runNumber << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(grpmagPath, runNumber);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for run " << runNumber;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      mVtx = ccdb->getForRun<o2::dataformats::MeanVertexObject>(mVtxPath, runNumber);
      LOG(info) << "Retrieved GRP for run " << runNumber << " with magnetic field of " << d_bz << " kZG";
    }

    // if TOF Nsigma desired
    if (doNSigmas) {
      nSigmaCalibObjects = ccdb->getForRun<TList>(nSigmaPath, runNumber);
      if (nSigmaCalibObjects) {
        LOGF(info, "loaded TList with this many objects: %i", nSigmaCalibObjects->GetEntries());

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
  void processV0Candidate(TCollision const& collision, TV0 const& v0, TTrack const& pTra, TTrack const& nTra)
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

    float deltaTimePositiveLambdaPi = -1e+6;
    float deltaTimeNegativeLambdaPi = -1e+6;
    float deltaTimePositiveLambdaPr = -1e+6;
    float deltaTimeNegativeLambdaPr = -1e+6;
    float deltaTimePositiveK0ShortPi = -1e+6;
    float deltaTimeNegativeK0ShortPi = -1e+6;

    float nSigmaPositiveLambdaPi = -1e+3;
    float nSigmaPositiveLambdaPr = -1e+3;
    float nSigmaNegativeLambdaPi = -1e+3;
    float nSigmaNegativeLambdaPr = -1e+3;
    float nSigmaPositiveK0ShortPi = -1e+3;
    float nSigmaNegativeK0ShortPi = -1e+3;

    float velocityPositivePr = velocity(posTrack.getP(), o2::constants::physics::MassProton);
    float velocityPositivePi = velocity(posTrack.getP(), o2::constants::physics::MassPionCharged);
    float velocityNegativePr = velocity(negTrack.getP(), o2::constants::physics::MassProton);
    float velocityNegativePi = velocity(negTrack.getP(), o2::constants::physics::MassPionCharged);

    float lengthPositive = findInterceptLength(posTrack, d_bz); // FIXME: tofPosition ok? adjust?
    float lengthNegative = findInterceptLength(negTrack, d_bz); // FIXME: tofPosition ok? adjust?
    float timePositivePr = lengthPositive / velocityPositivePr;
    float timePositivePi = lengthPositive / velocityPositivePi;
    float timeNegativePr = lengthNegative / velocityNegativePr;
    float timeNegativePi = lengthNegative / velocityNegativePi;

    if (pTra.hasTOF() && lengthPositive > 0) {
      deltaTimePositiveLambdaPr = (pTra.tofSignal() - pTra.tofEvTime()) - (timeLambda + timePositivePr);
      deltaTimePositiveLambdaPi = (pTra.tofSignal() - pTra.tofEvTime()) - (timeLambda + timePositivePi);
      deltaTimePositiveK0ShortPi = (pTra.tofSignal() - pTra.tofEvTime()) - (timeK0Short + timePositivePi);
    }
    if (nTra.hasTOF() && lengthNegative > 0) {
      deltaTimeNegativeLambdaPr = (nTra.tofSignal() - nTra.tofEvTime()) - (timeLambda + timeNegativePr);
      deltaTimeNegativeLambdaPi = (nTra.tofSignal() - nTra.tofEvTime()) - (timeLambda + timeNegativePi);
      deltaTimeNegativeK0ShortPi = (nTra.tofSignal() - nTra.tofEvTime()) - (timeK0Short + timeNegativePi);
    }

    if (doQA) {
      // calculate and pack properties for QA purposes
      int posProperties = 0;
      if (lengthPositive > 0)
        posProperties = posProperties | (static_cast<int>(1) << kLength);
      if (pTra.hasTOF())
        posProperties = posProperties | (static_cast<int>(1) << kHasTOF);
      int negProperties = 0;
      if (lengthNegative > 0)
        negProperties = negProperties | (static_cast<int>(1) << kLength);
      if (nTra.hasTOF())
        negProperties = negProperties | (static_cast<int>(1) << kHasTOF);

      histos.fill(HIST("h2dPositiveTOFProperties"), v0.pt(), posProperties);
      histos.fill(HIST("h2dNegativeTOFProperties"), v0.pt(), negProperties);
    }

    float deltaDecayTimeLambda = -10e+4;
    float deltaDecayTimeAntiLambda = -10e+4;
    float deltaDecayTimeK0Short = -10e+4;
    if (nTra.hasTOF() && pTra.hasTOF() > 0 && lengthPositive > 0 && lengthNegative > 0) { // does not depend on event time
      deltaDecayTimeLambda = (pTra.tofSignal() - timePositivePr) - (nTra.tofSignal() - timeNegativePi);
      deltaDecayTimeAntiLambda = (pTra.tofSignal() - timePositivePi) - (nTra.tofSignal() - timeNegativePr);
      deltaDecayTimeK0Short = (pTra.tofSignal() - timePositivePi) - (nTra.tofSignal() - timeNegativePi);
    }

    // calculate betas

    float evTimeMean = 0.5f * (pTra.tofEvTime() + nTra.tofEvTime());
    float decayTimeLambda = 0.5f * ((pTra.tofSignal() - timePositivePr) + (nTra.tofSignal() - timeNegativePi)) - evTimeMean;
    float decayTimeAntiLambda = 0.5f * ((pTra.tofSignal() - timePositivePi) + (nTra.tofSignal() - timeNegativePr)) - evTimeMean;
    float decayTimeK0Short = 0.5f * ((pTra.tofSignal() - timePositivePi) + (nTra.tofSignal() - timeNegativePi)) - evTimeMean;

    float betaLambda = -1e+6;
    float betaAntiLambda = -1e+6;
    float betaK0Short = -1e+6;

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
    if (doNSigmas) {
      // sweep through all viable hypotheses and produce N-sigma

      if (deltaTimePositiveLambdaPi > -1e+5)
        nSigmaPositiveLambdaPi = (deltaTimePositiveLambdaPi - hMeanPosLaPi->Interpolate(v0.pt())) / hSigmaPosLaPi->Interpolate(v0.pt());
      if (deltaTimePositiveLambdaPr > -1e+5)
        nSigmaPositiveLambdaPr = (deltaTimePositiveLambdaPr - hMeanPosLaPr->Interpolate(v0.pt())) / hSigmaPosLaPr->Interpolate(v0.pt());
      if (deltaTimeNegativeLambdaPi > -1e+5)
        nSigmaNegativeLambdaPi = (deltaTimeNegativeLambdaPi - hMeanNegLaPi->Interpolate(v0.pt())) / hSigmaNegLaPi->Interpolate(v0.pt());
      if (deltaTimeNegativeLambdaPr > -1e+5)
        nSigmaNegativeLambdaPr = (deltaTimeNegativeLambdaPr - hMeanNegLaPr->Interpolate(v0.pt())) / hSigmaNegLaPr->Interpolate(v0.pt());
      if (deltaTimePositiveK0ShortPi > -1e+5)
        nSigmaPositiveK0ShortPi = (deltaTimePositiveK0ShortPi - hMeanPosK0Pi->Interpolate(v0.pt())) / hSigmaPosK0Pi->Interpolate(v0.pt());
      if (deltaTimeNegativeK0ShortPi > -1e+5)
        nSigmaNegativeK0ShortPi = (deltaTimeNegativeK0ShortPi - hMeanNegK0Pi->Interpolate(v0.pt())) / hSigmaNegK0Pi->Interpolate(v0.pt());

      v0tofnsigmas(
        nSigmaPositiveLambdaPr, nSigmaNegativeLambdaPi,
        nSigmaNegativeLambdaPr, nSigmaPositiveLambdaPi,
        nSigmaPositiveK0ShortPi, nSigmaNegativeK0ShortPi);
    }

    if (doQA) {
      if (pTra.hasTOF()) {
        histos.fill(HIST("h2dProtonMeasuredVsExpected"),
                    (timeLambda + timePositivePr),
                    (pTra.tofSignal() - pTra.tofEvTime()));
        if (v0.v0cosPA() > qaCosPA && v0.dcaV0daughters() < qaDCADau) {
          if (std::abs(v0.mLambda() - 1.115683) < qaMassWindow && fabs(pTra.tpcNSigmaPr()) < qaTPCNSigma && fabs(nTra.tpcNSigmaPi()) < qaTPCNSigma) {
            histos.fill(HIST("h2dDeltaTimePositiveLambdaPr"), v0.pt(), v0.eta(), deltaTimePositiveLambdaPr);
            if (doQANSigma)
              histos.fill(HIST("h2dNSigmaPositiveLambdaPr"), v0.pt(), nSigmaPositiveLambdaPr);
          }
          if (std::abs(v0.mAntiLambda() - 1.115683) < qaMassWindow && fabs(pTra.tpcNSigmaPi()) < qaTPCNSigma && fabs(nTra.tpcNSigmaPr()) < qaTPCNSigma) {
            histos.fill(HIST("h2dDeltaTimePositiveLambdaPi"), v0.pt(), v0.eta(), deltaTimePositiveLambdaPi);
            if (doQANSigma)
              histos.fill(HIST("h2dNSigmaPositiveLambdaPi"), v0.pt(), nSigmaPositiveLambdaPi);
          }
          if (std::abs(v0.mK0Short() - 0.497) < qaMassWindow && fabs(pTra.tpcNSigmaPi()) < qaTPCNSigma && fabs(nTra.tpcNSigmaPi()) < qaTPCNSigma) {
            histos.fill(HIST("h2dDeltaTimePositiveK0ShortPi"), v0.pt(), v0.eta(), deltaTimePositiveK0ShortPi);
            if (doQANSigma)
              histos.fill(HIST("h2dNSigmaPositiveK0ShortPi"), v0.pt(), nSigmaPositiveK0ShortPi);
          }
        }
      }

      if (nTra.hasTOF()) {
        histos.fill(HIST("h2dPionMeasuredVsExpected"),
                    (timeLambda + timeNegativePi),
                    (nTra.tofSignal() - nTra.tofEvTime()));
        if (v0.v0cosPA() > qaCosPA && v0.dcaV0daughters() < qaDCADau) {
          if (std::abs(v0.mLambda() - 1.115683) < qaMassWindow && fabs(pTra.tpcNSigmaPr()) < qaTPCNSigma && fabs(nTra.tpcNSigmaPi()) < qaTPCNSigma) {
            histos.fill(HIST("h2dDeltaTimeNegativeLambdaPi"), v0.pt(), v0.eta(), deltaTimeNegativeLambdaPi);
            if (doQANSigma)
              histos.fill(HIST("h2dNSigmaNegativeLambdaPi"), v0.pt(), nSigmaNegativeLambdaPi);
          }
          if (std::abs(v0.mAntiLambda() - 1.115683) < qaMassWindow && fabs(pTra.tpcNSigmaPi()) < qaTPCNSigma && fabs(nTra.tpcNSigmaPr()) < qaTPCNSigma) {
            histos.fill(HIST("h2dDeltaTimeNegativeLambdaPr"), v0.pt(), v0.eta(), deltaTimeNegativeLambdaPr);
            if (doQANSigma)
              histos.fill(HIST("h2dNSigmaNegativeLambdaPr"), v0.pt(), nSigmaNegativeLambdaPr);
          }
          if (std::abs(v0.mK0Short() - 0.497) < qaMassWindow && fabs(pTra.tpcNSigmaPi()) < qaTPCNSigma && fabs(nTra.tpcNSigmaPi()) < qaTPCNSigma) {
            histos.fill(HIST("h2dDeltaTimeNegativeK0ShortPi"), v0.pt(), v0.eta(), deltaTimeNegativeK0ShortPi);
            if (doQANSigma)
              histos.fill(HIST("h2dNSigmaNegativeK0ShortPi"), v0.pt(), nSigmaNegativeK0ShortPi);
          }
        }
      }
      // delta lambda decay time
      histos.fill(HIST("h2dLambdaDeltaDecayTime"), v0.pt(), deltaDecayTimeLambda);
    }
  }

  void processStandardData(aod::Collisions const& collisions, V0OriginalDatas const& V0s, TracksWithAllExtras const&, aod::BCsWithTimestamps const& /*bcs*/)
  {
    // Fire up CCDB with first collision in record. If no collisions, bypass
    if (useCustomRunNumber || collisions.size() < 1) {
      initCCDB(manualRunNumber);
    } else {
      auto collision = collisions.begin();
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc.runNumber());
    }

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
      processV0Candidate(primaryVertex, V0, pTra, nTra);
    }
  }

  void processDerivedData(soa::Join<aod::StraCollisions, aod::StraStamps> const& collisions, V0DerivedDatas const& V0s, dauTracks const&)
  {
    // Fire up CCDB with first collision in record. If no collisions, bypass
    if (useCustomRunNumber || collisions.size() < 1) {
      initCCDB(manualRunNumber);
    } else {
      auto collision = collisions.begin();
      initCCDB(collision.runNumber());
    }

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
      processV0Candidate(primaryVertex, V0, pTra, nTra);
    }
  }

  PROCESS_SWITCH(lambdakzeropid, processStandardData, "Process standard data", true);
  PROCESS_SWITCH(lambdakzeropid, processDerivedData, "Process derived data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzeropid>(cfgc)};
}
