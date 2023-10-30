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
// O2 includes

/// \file HFFilterPrepareMLSamples.cxx
/// \brief task for trainings of ML models to be used in the HFFilter.cxx task
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Marcel Lesch <marcel.lesch@tum.de>, TUM
/// \author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
/// \author Biao Zhang <biao.zhang@cern.ch>, CCNU

#include <onnxruntime/core/session/experimental_onnxruntime_cxx_api.h> // needed for HFFilterHelpers, to be fixed

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/trackUtilities.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "EventFiltering/PWGHF/HFFilterHelpers.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod::hffilters;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfFilterPrepareMlSamples { // Main struct

  Produces<aod::HFTrigTrain2P> train2P;
  Produces<aod::HFTrigTrain3P> train3P;

  // parameters for production of training samples
  Configurable<bool> fillSignal{"fillSignal", true, "Flag to fill derived tables with signal for ML trainings"};
  Configurable<bool> fillBackground{"fillBackground", true, "Flag to fill derived tables with background for ML trainings"};
  Configurable<float> downSampleBkgFactor{"downSampleBkgFactor", 1., "Fraction of background candidates to keep for ML trainings"};

  // CCDB configuration
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  // material correction for track propagation
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
  int currentRun = 0; // needed to detect if the run changed and trigger update of calibrations etc.

  void init(InitContext&)
  {
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(url);
  }

  using BigTracksMCPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::McTrackLabels>;

  void process(aod::Hf2Prongs const& cand2Prongs,
               aod::Hf3Prongs const& cand3Prongs,
               aod::McParticles const& mcParticles,
               soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions,
               BigTracksMCPID const&,
               aod::BCsWithTimestamps const&)
  {
    for (const auto& cand2Prong : cand2Prongs) { // start loop over 2 prongs

      auto thisCollId = cand2Prong.collisionId();
      auto collision = collisions.rawIteratorAt(thisCollId);
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();

      if (currentRun != bc.runNumber()) {
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag, bc.timestamp());
        o2::base::Propagator::initFieldFromGRP(grpo);
        currentRun = bc.runNumber();
      }

      auto trackPos = cand2Prong.prong0_as<BigTracksMCPID>(); // positive daughter
      auto trackNeg = cand2Prong.prong1_as<BigTracksMCPID>(); // negative daughter

      auto trackParPos = getTrackPar(trackPos);
      auto trackParNeg = getTrackPar(trackNeg);
      o2::gpu::gpustd::array<float, 2> dcaPos{trackPos.dcaXY(), trackPos.dcaZ()};
      o2::gpu::gpustd::array<float, 2> dcaNeg{trackNeg.dcaXY(), trackNeg.dcaZ()};
      std::array<float, 3> pVecPos{trackPos.px(), trackPos.py(), trackPos.pz()};
      std::array<float, 3> pVecNeg{trackNeg.px(), trackNeg.py(), trackNeg.pz()};
      if (trackPos.collisionId() != thisCollId) {
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParPos, 2.f, noMatCorr, &dcaPos);
        getPxPyPz(trackParPos, pVecPos);
      }
      if (trackNeg.collisionId() != thisCollId) {
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParNeg, 2.f, noMatCorr, &dcaNeg);
        getPxPyPz(trackParNeg, pVecNeg);
      }

      auto pVec2Prong = RecoDecay::pVec(pVecPos, pVecNeg);
      auto pt2Prong = RecoDecay::pt(pVec2Prong);

      auto invMassD0 = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massPi, massKa});
      auto invMassD0bar = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massKa, massPi});

      int8_t sign = 0;
      int8_t flag = RecoDecay::OriginType::None;

      // D0(bar) → π± K∓
      bool isInCorrectColl{false};
      auto indexRec = RecoDecay::getMatchedMCRec(mcParticles, std::array{trackPos, trackNeg}, pdg::Code::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign);
      if (indexRec > -1) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        flag = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
        isInCorrectColl = (collision.mcCollisionId() == particle.mcCollisionId());
        if (flag < RecoDecay::OriginType::Prompt) {
          continue;
        }
      }

      float pseudoRndm = trackPos.pt() * 1000. - (int64_t)(trackPos.pt() * 1000);
      if ((fillSignal && indexRec > -1) || (fillBackground && indexRec < 0 && pseudoRndm < downSampleBkgFactor)) {
        train2P(invMassD0, invMassD0bar, pt2Prong, trackParPos.getPt(), dcaPos[0], dcaPos[1], trackPos.tpcNSigmaPi(), trackPos.tpcNSigmaKa(), trackPos.tofNSigmaPi(), trackPos.tofNSigmaKa(),
                trackParNeg.getPt(), dcaNeg[0], dcaNeg[1], trackNeg.tpcNSigmaPi(), trackNeg.tpcNSigmaKa(), trackNeg.tofNSigmaPi(), trackNeg.tofNSigmaKa(), flag, isInCorrectColl);
      }
    } // end loop over 2-prong candidates

    for (const auto& cand3Prong : cand3Prongs) { // start loop over 3 prongs

      auto thisCollId = cand3Prong.collisionId();
      auto collision = collisions.rawIteratorAt(thisCollId);
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();

      if (currentRun != bc.runNumber()) {
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag, bc.timestamp());
        o2::base::Propagator::initFieldFromGRP(grpo);
        currentRun = bc.runNumber();
      }

      auto trackFirst = cand3Prong.prong0_as<BigTracksMCPID>();  // first daughter
      auto trackSecond = cand3Prong.prong1_as<BigTracksMCPID>(); // second daughter
      auto trackThird = cand3Prong.prong2_as<BigTracksMCPID>();  // third daughter
      auto arrayDaughters = std::array{trackFirst, trackSecond, trackThird};

      auto trackParFirst = getTrackPar(trackFirst);
      auto trackParSecond = getTrackPar(trackSecond);
      auto trackParThird = getTrackPar(trackThird);
      o2::gpu::gpustd::array<float, 2> dcaFirst{trackFirst.dcaXY(), trackFirst.dcaZ()};
      o2::gpu::gpustd::array<float, 2> dcaSecond{trackSecond.dcaXY(), trackSecond.dcaZ()};
      o2::gpu::gpustd::array<float, 2> dcaThird{trackThird.dcaXY(), trackThird.dcaZ()};
      std::array<float, 3> pVecFirst{trackFirst.px(), trackFirst.py(), trackFirst.pz()};
      std::array<float, 3> pVecSecond{trackSecond.px(), trackSecond.py(), trackSecond.pz()};
      std::array<float, 3> pVecThird{trackThird.px(), trackThird.py(), trackThird.pz()};
      if (trackFirst.collisionId() != thisCollId) {
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParFirst, 2.f, noMatCorr, &dcaFirst);
        getPxPyPz(trackParFirst, pVecFirst);
      }
      if (trackSecond.collisionId() != thisCollId) {
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParSecond, 2.f, noMatCorr, &dcaSecond);
        getPxPyPz(trackParSecond, pVecSecond);
      }
      if (trackThird.collisionId() != thisCollId) {
        o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackParThird, 2.f, noMatCorr, &dcaThird);
        getPxPyPz(trackParThird, pVecThird);
      }

      auto pVec3Prong = RecoDecay::pVec(pVecFirst, pVecSecond, pVecThird);
      auto pt3Prong = RecoDecay::pt(pVec3Prong);

      auto invMassDplus = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massPi});

      auto invMassDsToKKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massKa, massKa, massPi});
      auto invMassDsToPiKK = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massKa});

      auto invMassLcToPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massProton, massKa, massPi});
      auto invMassLcToPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massProton});

      auto invMassXicToPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massProton, massKa, massPi});
      auto invMassXicToPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massKa, massProton});

      float deltaMassKKFirst = -1.f;
      float deltaMassKKSecond = -1.f;
      if (TESTBIT(cand3Prong.hfflag(), o2::aod::hf_cand_3prong::DecayType::DsToKKPi)) {
        deltaMassKKFirst = std::abs(RecoDecay::m(std::array{pVecFirst, pVecSecond}, std::array{massKa, massKa}) - massPhi);
        deltaMassKKSecond = std::abs(RecoDecay::m(std::array{pVecThird, pVecSecond}, std::array{massKa, massKa}) - massPhi);
      }
      int8_t sign = 0;
      int8_t flag = RecoDecay::OriginType::None;
      int8_t channel = -1;

      // D± → π± K∓ π±
      auto indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdg::Code::kDPlus, std::array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
      if (indexRec >= 0) {
        channel = kDplus;
      }
      if (indexRec < 0) {
        // Ds± → K± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdg::Code::kDS, std::array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec >= 0) {
          channel = kDs;
        }
      }
      if (indexRec < 0) {
        // Λc± → p± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdg::Code::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec >= 0) {
          channel = kLc;
        }
      }
      if (indexRec < 0) {
        // Ξc± → p± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(mcParticles, arrayDaughters, pdg::Code::kXiCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec >= 0) {
          channel = kXic;
        }
      }

      bool isInCorrectColl{false};
      if (indexRec > -1) {
        auto particle = mcParticles.rawIteratorAt(indexRec);
        flag = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
        isInCorrectColl = (collision.mcCollisionId() == particle.mcCollisionId());
        if (flag < RecoDecay::OriginType::Prompt) {
          continue;
        }
      }

      float pseudoRndm = trackFirst.pt() * 1000. - (int64_t)(trackFirst.pt() * 1000);
      if ((fillSignal && indexRec > -1) || (fillBackground && indexRec < 0 && pseudoRndm < downSampleBkgFactor)) {
        train3P(invMassDplus, invMassDsToKKPi, invMassDsToPiKK, invMassLcToPKPi, invMassLcToPiKP, invMassXicToPKPi, invMassXicToPiKP, pt3Prong, deltaMassKKFirst, deltaMassKKSecond,
                trackParFirst.getPt(), dcaFirst[0], dcaFirst[1], trackFirst.tpcNSigmaPi(), trackFirst.tpcNSigmaKa(), trackFirst.tpcNSigmaPr(), trackFirst.tofNSigmaPi(), trackFirst.tofNSigmaKa(), trackFirst.tofNSigmaPr(),
                trackParSecond.getPt(), dcaSecond[0], dcaSecond[1], trackSecond.tpcNSigmaPi(), trackSecond.tpcNSigmaKa(), trackSecond.tpcNSigmaPr(), trackSecond.tofNSigmaPi(), trackSecond.tofNSigmaKa(), trackSecond.tofNSigmaPr(),
                trackParThird.getPt(), dcaThird[0], dcaThird[1], trackThird.tpcNSigmaPi(), trackThird.tpcNSigmaKa(), trackThird.tpcNSigmaPr(), trackThird.tofNSigmaPi(), trackThird.tofNSigmaKa(), trackThird.tofNSigmaPr(),
                flag, channel, cand3Prong.hfflag(), isInCorrectColl);
      }
    } // end loop over 3-prong candidates
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{

  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfFilterPrepareMlSamples>(cfg));

  return workflow;
}
