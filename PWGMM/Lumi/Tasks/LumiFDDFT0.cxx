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
// author: akhuntia@cern.ch

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "CommonUtils/NameConf.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"

#include "DetectorsVertexing/PVertexer.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/Vertex.h"

#include "DataFormatsFDD/Digit.h"
#include "DataFormatsFIT/Triggers.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "CommonConstants/GeomConstants.h"
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "DataFormatsCalibration/MeanVertexObject.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using BCsWithTimestamps = soa::Join<aod::BCs, aod::Timestamps>;
int nBCsPerOrbit = 3564;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, uint64_t);
DECLARE_SOA_COLUMN(VertexX, vertexX, double);
DECLARE_SOA_COLUMN(VertexY, vertexY, double);
DECLARE_SOA_COLUMN(VertexZ, vertexZ, double);
DECLARE_SOA_COLUMN(VertexXX, vertexXX, double);
DECLARE_SOA_COLUMN(VertexYY, vertexYY, double);
DECLARE_SOA_COLUMN(VertexXY, vertexXY, double);
DECLARE_SOA_COLUMN(GlobalBC, globalBC, uint64_t);
DECLARE_SOA_COLUMN(VertexChi2, vertexChi2, double);
DECLARE_SOA_COLUMN(NContrib, nContrib, int);

// Information for FDD
DECLARE_SOA_COLUMN(isFDD, isfdd, bool);
DECLARE_SOA_COLUMN(TCMTriggerFDD, tcmTriggerfdd, uint8_t);
DECLARE_SOA_COLUMN(TimeAFDD, timeAfdd, double);
DECLARE_SOA_COLUMN(TimeCFDD, timeCfdd, double);
DECLARE_SOA_COLUMN(ChargeAFDD, chargeAfdd, double);
DECLARE_SOA_COLUMN(ChargeCFDD, chargeCfdd, double);

// Information for FT0
DECLARE_SOA_COLUMN(isFT0, isft0, bool);
DECLARE_SOA_COLUMN(TCMTriggerFT0, tcmTriggerft0, uint8_t);
DECLARE_SOA_COLUMN(TimeAFT0, timeAft0, double);
DECLARE_SOA_COLUMN(TimeCFT0, timeCft0, double);
DECLARE_SOA_COLUMN(ChargeAFT0, chargeAft0, double);
DECLARE_SOA_COLUMN(ChargeCFT0, chargeCft0, double);

} // namespace full
DECLARE_SOA_TABLE(EventInfo, "AOD", "EventInfo", full::TimeStamp, full::VertexX,
                  full::VertexY, full::VertexZ, full::GlobalBC,
                  full::VertexChi2, full::NContrib,
                  full::isFDD, full::TCMTriggerFDD,
                  full::TimeAFDD, full::TimeCFDD,
                  full::ChargeAFDD, full::ChargeCFDD,
                  full::isFT0, full::TCMTriggerFT0,
                  full::TimeAFT0, full::TimeCFT0,
                  full::ChargeAFT0, full::ChargeCFT0);

DECLARE_SOA_TABLE(EventInfoFDD, "AOD", "EventInfoFDD",
                  full::TimeStamp, full::GlobalBC,
                  full::TCMTriggerFDD, full::TimeAFDD,
                  full::TimeCFDD, full::ChargeAFDD,
                  full::ChargeCFDD);

DECLARE_SOA_TABLE(EventInfoFT0, "AOD", "EventInfoFT0",
                  full::TimeStamp, full::GlobalBC,
                  full::TCMTriggerFT0, full::TimeAFT0,
                  full::TimeCFT0, full::ChargeAFT0,
                  full::ChargeCFT0);

} // namespace o2::aod

struct LumiFDDFT0 {
  Produces<o2::aod::EventInfo> rowEventInfo;
  Produces<o2::aod::EventInfoFDD> rowEventInfofdd;
  Produces<o2::aod::EventInfoFT0> rowEventInfoft0;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  const char* ccdbpath_grp = "GLO/Config/GRPMagField";
  const char* ccdburl = "http://alice-ccdb.cern.ch";
  int mRunNumber;

  Configurable<uint64_t> fttimestamp{"fttimestamp", 1668080173000, "First time of time stamp"};
  Configurable<int> nContribMax{"nContribMax", 2500, "Maximum number of contributors"};
  Configurable<int> nContribMin{"nContribMin", 10, "Minimum number of contributors"};

  HistogramRegistry histos{
    "histos",
    {
      {"vertexx", "", {HistType::kTH1F, {{1000, -1, 1, "x"}}}},                                         //
      {"vertexy", "", {HistType::kTH1F, {{1000, -1, 1, "y"}}}},                                         //
      {"timestamp", "", {HistType::kTH1F, {{20000, 0, 2e7, "t"}}}},                                     //
      {"vertexx_timestamp", "", {HistType::kTH2F, {{10000, 130e9, 140e9, "t"}, {2000, -1, 1, "x"}}}},   //
      {"vertexy_timestamp", "", {HistType::kTH2F, {{10000, 130e9, 140e9, "t"}, {2000, -1, 1, "y"}}}},   //
      {"chisquare", "", {HistType::kTH1F, {{1000, 0, 100, "#chi^{2}"}}}},                               //
      {"vertexx_Refitted", "", {HistType::kTH1F, {{1000, -1, 1, "x"}}}},                                //
      {"vertexy_Refitted", "", {HistType::kTH1F, {{1000, -1, 1, "y"}}}},                                //
      {"vertexx_Refitted_timestamp", "", {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {2000, -1, 1, "x"}}}}, //
      {"vertexy_Refitted_timestamp", "", {HistType::kTH2F, {{1000, 0, 4e6, "t"}, {2000, -1, 1, "y"}}}}, //
      {"chisquare_Refitted", "", {HistType::kTH1F, {{1000, 0, 100, "#chi^{2}"}}}},                      //
      {"vertexx_Refitted_vertexx", "", {HistType::kTH2F, {{1000, -1, 1, "x"}, {1000, -1, 1, "rx"}}}},   //
      {"vertexy_Refitted_vertexy", "", {HistType::kTH2F, {{1000, -1, 1, "y"}, {1000, -1, 1, "ry"}}}}    //
    }};

  HistogramRegistry histoslite{
    "histoslite",
    {
      {"BCFDD", "", {HistType::kTH1F, {{nBCsPerOrbit + 1, -0.5f, nBCsPerOrbit + 0.5f, "x"}}}}, //
      {"BCFT0", "", {HistType::kTH1F, {{nBCsPerOrbit + 1, -0.5f, nBCsPerOrbit + 0.5f, "x"}}}}, //
    }};

  bool doPVrefit = true;

  void init(InitContext&)
  {
    if (doprocessLite == true && doprocessFull == true) {
      LOG(fatal) << "Select one process function from  processLite and processFull";
    }

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
    mRunNumber = 0;
  }

  void processFull(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::FDDs const& fdds, aod::FT0s const& ft0s, aod::BCsWithTimestamps const&,
                   o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov,
                                 o2::aod::TracksExtra> const& unfiltered_tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    Long64_t relTS = bc.timestamp() - fttimestamp;
    Long64_t globalBC = bc.globalBC();
    std::vector<int64_t> vec_globID_contr = {};
    std::vector<o2::track::TrackParCov> vec_TrkContributos = {};

    int nContrib = 0;
    int nNonContrib = 0;

    for (const auto& unfiltered_track : unfiltered_tracks) {
      if (!unfiltered_track.hasITS()) {
        nNonContrib++;
        continue;
      }
      if (unfiltered_track.pt() < 0.8 || unfiltered_track.itsNCls() < 5) {
        nNonContrib++;
        continue;
      }
      vec_globID_contr.push_back(unfiltered_track.globalIndex());
      vec_TrkContributos.push_back(getTrackParCov(unfiltered_track));
      nContrib++;
    }
    std::vector<bool> vec_useTrk_PVrefit(vec_globID_contr.size(), true);

    if (mRunNumber != bc.runNumber()) {
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbpath_grp, bc.timestamp());
      if (grpo != nullptr) {
        o2::base::Propagator::initFieldFromGRP(grpo);
      } else {
        LOGF(fatal,
             "GRP object is not available in CCDB for run=%d at timestamp=%llu",
             bc.runNumber(), bc.timestamp());
      }
      mRunNumber = bc.runNumber();
    }

    o2::dataformats::VertexBase Pvtx;
    Pvtx.setX(collision.posX());
    Pvtx.setY(collision.posY());
    Pvtx.setZ(collision.posZ());
    Pvtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(),
                collision.covXZ(), collision.covYZ(), collision.covZZ());
    // configure PVertexer
    o2::dataformats::VertexBase PVbase_recalculated;
    o2::vertexing::PVertexer vertexer;
    o2::conf::ConfigurableParam::updateFromString(
      "pvertexer.useMeanVertexConstraint=false"); // we want to refit w/o
                                                  // MeanVertex constraint
    vertexer.init();
    bool PVrefit_doable = vertexer.prepareVertexRefit(vec_TrkContributos, Pvtx);
    double chi2 = -1.;
    double refitX = -999.;
    double refitY = -999.;
    double refitZ = -999.;
    double refitXX = -999.;
    double refitYY = -999.;
    double refitXY = -999.;

    double timeaFDD = -999.;
    double timecFDD = -999.;
    double chargeaFDD = 0.;
    double chargecFDD = 0.;
    uint8_t mTriggerFDD = 0;

    double timeaFT0 = -999.;
    double timecFT0 = -999.;
    double chargeaFT0 = 0.;
    double chargecFT0 = 0.;
    uint8_t mTriggerFT0 = 0;

    if (doPVrefit && PVrefit_doable) {
      auto Pvtx_refitted = vertexer.refitVertex(vec_useTrk_PVrefit, Pvtx);
      chi2 = Pvtx_refitted.getChi2();
      refitX = Pvtx_refitted.getX();
      refitY = Pvtx_refitted.getY();
      refitZ = Pvtx_refitted.getZ();
      refitXX = Pvtx_refitted.getSigmaX2();
      refitYY = Pvtx_refitted.getSigmaY2();
      refitXY = Pvtx_refitted.getSigmaXY();

      // now get information for FDD
      if (collision.has_foundFDD()) {
        auto fdd = collision.foundFDD();
        mTriggerFDD = fdd.triggerMask();
        timeaFDD = fdd.timeA();
        timecFDD = fdd.timeC();
        for (auto amplitude : fdd.chargeA()) {
          chargeaFDD += amplitude;
        }
        for (auto amplitude : fdd.chargeC()) {
          chargecFDD += amplitude;
        }
      } // fdd

      if (collision.has_foundFT0()) {
        auto ft0 = collision.foundFT0();
        mTriggerFT0 = ft0.triggerMask();
        timeaFT0 = ft0.timeA();
        timecFT0 = ft0.timeC();
        for (auto amplitude : ft0.amplitudeA()) {
          chargeaFT0 += amplitude;
        }

        for (auto amplitude : ft0.amplitudeC()) {
          chargecFT0 += amplitude;
        }
      } // ft0

    } // pv refit
    rowEventInfo(relTS, refitX, refitY, refitZ, globalBC, chi2, nContrib, collision.has_foundFDD(),
                 mTriggerFDD, timeaFDD, timecFDD, chargeaFDD, chargecFDD, collision.has_foundFT0(), mTriggerFT0, timeaFT0,
                 timecFT0, chargeaFT0, chargecFT0);

    histos.fill(HIST("chisquare_Refitted"), chi2);
    if (nContrib > nContribMin && nContrib < nContribMax &&
        (chi2 / nContrib) < 4.0 && chi2 > 0) {
      histos.fill(HIST("vertexx_Refitted"), refitX);
      histos.fill(HIST("vertexy_Refitted"), refitY);

      histos.fill(HIST("vertexx_Refitted_timestamp"), relTS, refitX);
      histos.fill(HIST("vertexy_Refitted_timestamp"), relTS, refitY);
    }
    histos.fill(HIST("chisquare"), collision.chi2());
    if (collision.chi2() / collision.numContrib() > 4)
      return;
    if (collision.numContrib() > nContribMax ||
        collision.numContrib() < nContribMin)
      return;
    histos.fill(HIST("vertexx"), collision.posX());
    histos.fill(HIST("vertexy"), collision.posY());
    histos.fill(HIST("timestamp"), relTS);

    histos.fill(HIST("vertexx_timestamp"), relTS, collision.posX());
    histos.fill(HIST("vertexy_timestamp"), relTS, collision.posY());

    if (nContrib > nContribMin && nContrib < nContribMax &&
        (chi2 / nContrib) < 4.0 && chi2 > 0) {
      histos.fill(HIST("vertexx_Refitted_vertexx"), collision.posX(), refitX);
      histos.fill(HIST("vertexy_Refitted_vertexy"), collision.posY(), refitY);
    }
    // need selections
  };
  PROCESS_SWITCH(LumiFDDFT0, processFull, "Process FDD", true);

  void processLite(aod::FDDs const& fdds, aod::FT0s const& ft0s, aod::BCsWithTimestamps const&)
  {

    // Scan over the FDD table and store charge and time along with globalBC
    for (auto& fdd : fdds) {
      auto bc = fdd.bc_as<BCsWithTimestamps>();
      if (!bc.timestamp())
        continue;
      if (mRunNumber != bc.runNumber()) {
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbpath_grp, bc.timestamp());
        if (grpo != nullptr) {
          o2::base::Propagator::initFieldFromGRP(grpo);
          std::cout << "run " << bc.runNumber() << std::endl;
        } else {
          LOGF(fatal,
               "GRP object is not available in CCDB for run=%d at timestamp=%llu",
               bc.runNumber(), bc.timestamp());
        }
        mRunNumber = bc.runNumber();
      }
      // get the relative ts and globalBC, which will help to sync FDD and FT0 data
      Long64_t relTS = bc.timestamp() - fttimestamp;
      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;
      histoslite.fill(HIST("BCFDD"), localBC);

      double chargeaFDD = 0.;
      double chargecFDD = 0.;

      for (auto amplitude : fdd.chargeA()) {
        chargeaFDD += amplitude;
      }
      for (auto amplitude : fdd.chargeC()) {
        chargecFDD += amplitude;
      }
      rowEventInfofdd(relTS, globalBC, fdd.triggerMask(), fdd.timeA(), fdd.timeC(), chargeaFDD, chargecFDD);
    } // end of fdd table

    // Scan over the FT0 table and store charge and time along with globalBC
    for (auto& ft0 : ft0s) {
      auto bc = ft0.bc_as<BCsWithTimestamps>();
      if (!bc.timestamp())
        continue;
      Long64_t relTS = bc.timestamp() - fttimestamp;
      Long64_t globalBC = bc.globalBC();
      int localBC = globalBC % nBCsPerOrbit;
      histoslite.fill(HIST("BCFT0"), localBC);
      double chargeaFT0 = 0.;
      double chargecFT0 = 0.;
      for (auto amplitude : ft0.amplitudeA()) {
        chargeaFT0 += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        chargecFT0 += amplitude;
      }
      rowEventInfoft0(relTS, globalBC, ft0.triggerMask(), ft0.timeA(), ft0.timeC(), chargeaFT0, chargecFT0);
    } // end of ft0 table
  };
  PROCESS_SWITCH(LumiFDDFT0, processLite, "Process FDD and FT0 info", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w{adaptAnalysisTask<LumiFDDFT0>(cfgc, TaskName{"LumiFDDFT0"})};
  return w;
}
