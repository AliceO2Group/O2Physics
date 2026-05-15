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

#include "tableHMPID.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <TGeoManager.h>

#include <HMPIDBase/Param.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_set>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct HmpidTableProducer {

  Produces<aod::HmpidAnalysis> hmpidAnalysis;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  const AxisSpec axisEvtCounter{1, 0, +1, ""};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  struct : ConfigurableGroup {
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the material LUT"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the aligned geometry"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of GRPMagField"};
  } ccdbConfig;

  // -----------------------------------------------------------------------
  // Quality configurables
  // -----------------------------------------------------------------------
  Configurable<bool> requireITS{"requireITS", true, "Require ITS track"};
  Configurable<bool> requireTPC{"requireTPC", true, "Require TPC track"};
  Configurable<bool> requireTOF{"requireTOF", true, "Require TOF track"};

  using CollisionCandidates = o2::soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFV0As>;

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra,
                                    aod::TracksDCA, aod::TrackSelection,
                                    aod::pidTPCFullPi, aod::pidTPCFullKa,
                                    aod::pidTPCFullPr, aod::pidTPCFullDe,
                                    aod::pidTOFFullPi, aod::pidTOFFullKa,
                                    aod::pidTOFFullPr, aod::pidTOFFullDe>;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdbConfig.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // event counters
    histos.add("eventCounter", "All events", kTH1F, {axisEvtCounter});
    histos.add("goodEventCounter", "Events passing sel8", kTH1F, {axisEvtCounter});
    histos.add("eventsHmpid", "Events with HMPID track", kTH1F, {axisEvtCounter});

    // chamber extraction quality checks - M1 propagation, M2 clusSize encoding, M3 hybrid configuration
    const AxisSpec axisCh{10, -1.5, 8.5, "Chamber"};
    histos.add("hChamberM1", "Chamber M1 (propagation)", kTH1F, {axisCh});
    histos.add("hChamberM2", "Chamber M2 (clusSize encoding)", kTH1F, {axisCh});
    histos.add("hChamberM3", "Chamber M3 (hybrid, in table)", kTH1F, {axisCh});

    histos.add("hChamberM1vsM2", "M1 vs M2; M2; M1", kTH2F, {axisCh, axisCh});
    histos.add("hChamberM3vsM2", "M3 vs M2; M2; M3", kTH2F, {axisCh, axisCh});

    histos.add("hClusSize", "Raw hmpidClusSize", kTH1F, {{500, -1.1e6, 1e6, "clusSize"}});
    histos.add("hClusSizeCorrupt", "Corrupt entries (<=0)", kTH1F, {{200, -1.1e6, 1., "clusSize"}});

    histos.add("hChamberAssignment",
               "Chamber assignment outcome; category; counts",
               kTH1F, {{4, -0.5, 3.5, ""}});
  }

  // -----------------------------------------------------------------------
  // CCDB initialisation per run
  // -----------------------------------------------------------------------
  int mCCDBRunNumber = 0;

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mCCDBRunNumber == bc.runNumber()) {
      return;
    }
    mCCDBRunNumber = bc.runNumber();

    auto grpMag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(
      ccdbConfig.grpmagPath, bc.timestamp());
    if (!grpMag) {
      LOGF(fatal, "initCCDB: GRPMagField not found at %s",
           ccdbConfig.grpmagPath.value.c_str());
      return;
    }
    o2::base::Propagator::initFieldFromGRP(grpMag);

    auto lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(
      ccdb->get<o2::base::MatLayerCylSet>(ccdbConfig.lutPath));
    o2::base::Propagator::Instance()->setMatLUT(lut);

    // Geometry needed by o2::hmpid::Param for mars2Lors
    if (!gGeoManager) {
      ccdb->get<TGeoManager>(ccdbConfig.geoPath);
    }

    LOGF(info, "initCCDB: run %d initialised (mag + LUT + geo)", bc.runNumber());
  }

  void propagateHelix(
    double len,
    double bz,
    int charge,
    std::array<double, 3>& x,
    std::array<double, 3>& p)
  {
    // Extrapolate track along simple helix in magnetic field

    constexpr double KB2C = 0.299792458e-3;

    double px = p[0];
    double py = p[1];
    double pz = p[2];

    double pmod = std::hypot(px, py, pz);
    double pt = std::hypot(px, py);

    if (pt < o2::constants::math::Almost0 || std::abs(bz) < o2::constants::math::Almost0) // straight-line tracks
    {
      x[0] += px / pmod * len;
      x[1] += py / pmod * len;
      x[2] += pz / pmod * len;

      return;
    }

    double a = -KB2C * bz * charge;
    double rho = a / pmod;

    double sinr = std::sin(rho * len);
    double cosr = std::cos(rho * len);

    x[0] += px * sinr / a - py * (1. - cosr) / a;
    x[1] += py * sinr / a + px * (1. - cosr) / a;
    x[2] += pz * len / pmod;

    double px0 = px;

    p[0] = px0 * cosr - py * sinr;
    p[1] = py * cosr + px0 * sinr;
    p[2] = pz;
  }

  bool intersectHelixPlane(
    double bz,
    int charge,
    const std::array<double, 3>& x0,
    const std::array<double, 3>& p0,
    const std::array<double, 3>& planePoint,
    const std::array<double, 3>& planeNormal,
    std::array<double, 3>& xOut,
    std::array<double, 3>& pOut)
  {
    // Intersect a helix with a plane defined by normals and points, using an iterative approach.

    double s =
      (planePoint[0] - x0[0]) * planeNormal[0] +
      (planePoint[1] - x0[1]) * planeNormal[1] +
      (planePoint[2] - x0[2]) * planeNormal[2];

    double dist = 99999., distPrev = dist;

    const double kConvThr = 0.00001;
    const int kMaxIter = 100;

    std::array<double, 3> x, p;

    for (int iter = 0; iter < kMaxIter; ++iter) {

      x = x0;
      p = p0;
      propagateHelix(s, bz, charge, x, p);

      dist =
        (x[0] - planePoint[0]) * planeNormal[0] +
        (x[1] - planePoint[1]) * planeNormal[1] +
        (x[2] - planePoint[2]) * planeNormal[2];

      if (std::abs(dist) >= std::abs(distPrev))
        return false;

      distPrev = dist;
      s -= dist;

      if (std::abs(dist) < kConvThr) {
        xOut = x;
        pOut = p;
        return true;
      }
    }

    return false;
  }

  int getHmpidChamber(
    std::array<double, 3> x,
    std::array<double, 3> p,
    double bz,
    int charge)
  {
    auto* param = o2::hmpid::Param::instance();

    for (int ch = o2::hmpid::Param::kMinCh; ch <= o2::hmpid::Param::kMaxCh; ++ch) {

      // Chamber module geometry: points and normals for radiators and PCs
      std::array<double, 3> pRad{}, pPc{};
      std::array<double, 3> nRad{}, nPc{};

      param->point(ch, pRad.data(), o2::hmpid::Param::kRad);
      param->norm(ch, nRad.data());
      param->point(ch, pPc.data(), o2::hmpid::Param::kPc);
      param->norm(ch, nPc.data());

      // Intersection track - radiator plane
      std::array<double, 3> xRad, pAtRad;

      if (!intersectHelixPlane(bz, charge, x, p, pRad, nRad, xRad, pAtRad))
        continue;

      // Intersection track - PC plane
      std::array<double, 3> xPc, pAtPc;

      if (!intersectHelixPlane(bz, charge, xRad, pAtRad, pPc, nPc, xPc, pAtPc))
        continue;

      double theta, phi;
      param->mars2LorsVec(ch, pAtRad.data(), theta, phi);

      double xL, yL;
      param->mars2Lors(ch, xPc.data(), xL, yL);

      // Use isInside to check Chamber intersected
      if (param->isInside(xL, yL, param->distCut()))
        return ch;
    }

    // No chamber intersected
    return -1;
  }

  void processEvent(CollisionCandidates::iterator const& col,
                    aod::BCsWithTimestamps const&)
  {
    histos.fill(HIST("eventCounter"), 0.5);
    if (col.sel8()) {
      histos.fill(HIST("goodEventCounter"), 0.5);
    }
    initCCDB(col.bc_as<aod::BCsWithTimestamps>());
  }
  PROCESS_SWITCH(HmpidTableProducer, processEvent, "Process event level", true);

  void processHmpid(
    aod::HMPIDs const& hmpids,
    TrackCandidates const&,
    CollisionCandidates const&,
    aod::BCsWithTimestamps const&)
  {
    static std::unordered_set<uint32_t> collisionsWithHmpid;

    for (auto const& t : hmpids) {

      const auto& globalTrack = t.track_as<TrackCandidates>();

      if (!globalTrack.has_collision())
        continue;

      const auto& col = globalTrack.collision_as<CollisionCandidates>();
      initCCDB(col.bc_as<aod::BCsWithTimestamps>());
      uint32_t collId = col.globalIndex();

      // Track quality selection
      if ((requireITS && !globalTrack.hasITS()) ||
          (requireTPC && !globalTrack.hasTPC()) ||
          (requireTOF && !globalTrack.hasTOF())) {
        continue;
      }

      if (collisionsWithHmpid.insert(collId).second) {
        histos.fill(HIST("eventsHmpid"), 0.5);
      }

      // clusSize diagnostics
      histos.fill(HIST("hClusSize"), t.hmpidClusSize());
      bool isCorrupt = (t.hmpidClusSize() <= 0);
      if (isCorrupt) {
        histos.fill(HIST("hClusSizeCorrupt"), t.hmpidClusSize());
      }

      // --- M2: clusSize encoding ---
      int chamberM2 = t.hmpidClusSize() / 1000000;
      histos.fill(HIST("hChamberM2"), chamberM2);

      // --- M1: propagation
      // obtain global coordinates
      double sinA = std::sin(globalTrack.alpha());
      double cosA = std::cos(globalTrack.alpha());

      std::array<double, 3> x = {
        globalTrack.x() * cosA - globalTrack.y() * sinA,
        globalTrack.x() * sinA + globalTrack.y() * cosA,
        static_cast<double>(globalTrack.z())};

      std::array<double, 3> p = {
        static_cast<double>(globalTrack.px()),
        static_cast<double>(globalTrack.py()),
        static_cast<double>(globalTrack.pz())};

      // int charge = (globalTrack.signed1Pt() > 0) ? +1 : -1;
      int16_t charge = globalTrack.sign();

      auto prop = o2::base::Propagator::Instance();
      double bz = static_cast<double>(prop->getNominalBz()); // positive sign

      int chamberM1 = getHmpidChamber(x, p, bz, charge);

      histos.fill(HIST("hChamberM1"), chamberM1);

      if (!isCorrupt) {
        histos.fill(HIST("hChamberM1vsM2"), chamberM2, chamberM1);
      }

      // --- M3: hybrid ---
      int chamberM3 = -1;
      if (!isCorrupt) {
        chamberM3 = chamberM2;
      } else {
        chamberM3 = chamberM1;
      }

      // Legend:
      // bin 0 = clusSize > 0,  chamber found   (M2 ok)
      // bin 1 = clusSize > 0,  chamber not found (non dovrebbe mai accadere)
      // bin 2 = clusSize <= 0, M1 recovery      (corrupt, M1 ok)
      // bin 3 = clusSize <= 0, M1 fails        (corrupt, skipped)

      if (!isCorrupt && chamberM3 >= 0) {
        histos.fill(HIST("hChamberAssignment"), 0.);
      } else if (!isCorrupt && chamberM3 < 0) {
        histos.fill(HIST("hChamberAssignment"), 1.);
      } else if (isCorrupt && chamberM3 >= 0) {
        histos.fill(HIST("hChamberAssignment"), 2.);
      } else {
        histos.fill(HIST("hChamberAssignment"), 3.);
      }

      histos.fill(HIST("hChamberM3"), chamberM3);
      histos.fill(HIST("hChamberM3vsM2"), chamberM2, chamberM3);

      // Skip track if chamber undetermined
      if (chamberM3 < 0) {
        continue;
      }

      // Fill photon charges
      float hmpidPhotsCharge2[o2::aod::kDimPhotonsCharge];
      for (int i = 0; i < o2::aod::kDimPhotonsCharge; i++) {
        hmpidPhotsCharge2[i] = t.hmpidPhotsCharge()[i];
      }

      // Fill output table
      hmpidAnalysis(
        t.hmpidSignal(), t.hmpidMom(),
        globalTrack.p(), t.hmpidXTrack(), t.hmpidYTrack(),
        t.hmpidXMip(), t.hmpidYMip(),
        t.hmpidNPhotons(), t.hmpidQMip(),
        (t.hmpidClusSize() % 1000000) / 1000,
        chamberM3,
        hmpidPhotsCharge2,
        globalTrack.eta(), globalTrack.phi(),
        globalTrack.px(), globalTrack.py(), globalTrack.pz(),
        globalTrack.itsNCls(), globalTrack.tpcNClsFound(), globalTrack.tpcNClsCrossedRows(),
        globalTrack.tpcChi2NCl(), globalTrack.itsChi2NCl(),
        globalTrack.dcaXY(), globalTrack.dcaZ(),
        globalTrack.tpcNSigmaPi(), globalTrack.tofNSigmaPi(),
        globalTrack.tpcNSigmaKa(), globalTrack.tofNSigmaKa(),
        globalTrack.tpcNSigmaPr(), globalTrack.tofNSigmaPr(),
        globalTrack.tpcNSigmaDe(), globalTrack.tofNSigmaDe(),
        col.centFV0A());

    } // end HMPID loop
  }
  PROCESS_SWITCH(HmpidTableProducer, processHmpid, "Process HMPID entries", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<HmpidTableProducer>(cfg)};
}
