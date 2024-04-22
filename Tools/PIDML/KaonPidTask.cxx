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

/// This task produces invariant mass vs. momentum and dEdX in TPC vs. momentum
/// for Kaons using ML PID from the PID ML ONNX Model.

#include <cmath>
#include <memory>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "TLorentzVector.h"
#include "TDatabasePDG.h"
#include "Framework/AnalysisDataModel.h"
#include "Tools/PIDML/pidOnnxModel.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"
#include "TMath.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
using MyCollisions = soa::Join<aod::Collisions,
                               aod::EvSels,
                               aod::Mults>;
using MyTracks = soa::Join<aod::FullTracks, aod::TracksExtra, aod::pidTOFbeta,
                           aod::TOFSignal, aod::TracksDCA>;
using MyCollision = MyCollisions::iterator;
using MyTrack = MyTracks::iterator;
} // namespace o2::aod

struct KaonPidTask {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = aod::track::collisionId;

  std::shared_ptr<PidONNXModel> pidModel; // creates a shared pointer to a new instance 'pidmodel'.
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> cfgZvtxCut{"cfgZvtxCut", 10, "Z vtx cut"};
  Configurable<float> cfgEtaCut{"cfgEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<float> cfgMaxPtCut{"cfgMaxPtCut", 3.0, "Max Pt cut"};
  Configurable<float> cfgMinPtCut{"cfgMinPtCut", 0.5, "Min Pt cut"};
  Configurable<float> cfgMinNSigmaTPCCut{"cfgMinNSigmaTPCCut", 3., "N-sigma TPC cut"};
  Configurable<float> cfgChargeCut{"cfgChargeCut", 0., "N-sigma TPC cut"};
  Configurable<std::string> cfgPathLocal{"local-path", ".", "base path to the local directory with ONNX models"};
  Configurable<std::string> cfgPathCCDB{"ccdb-path", "Users/m/mkabus/PIDML", "base path to the CCDB directory with ONNX models"};
  Configurable<std::string> cfgCCDBURL{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<int> cfgPid{"pid", 321, "PID to predict"};
  Configurable<double> cfgCertainty{"certainty", 0.5, "Minimum certainty above which the model accepts a particular type of particle"};
  Configurable<uint32_t> cfgDetector{"detector", kTPCTOFTRD, "What detectors to use: 0: TPC only, 1: TPC + TOF, 2: TPC + TOF + TRD"};
  Configurable<uint64_t> cfgTimestamp{"timestamp", 0, "Fixed timestamp"};
  Configurable<bool> cfgUseCCDB{"useCCDB", false, "Whether to autofetch ML model from CCDB. If false, local file will be used."};

  o2::ccdb::CcdbApi ccdbApi;

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgZvtxCut);
  Filter trackFilter = (nabs(aod::track::eta) < cfgEtaCut) && (aod::track::pt > cfgMinPtCut) && (aod::track::pt < cfgMaxPtCut);

  // Applying filters
  using MyFilteredCollisions = soa::Filtered<o2::aod::MyCollisions>;
  using MyFilteredCollision = MyFilteredCollisions::iterator;

  Partition<o2::aod::MyTracks> positive = (nabs(aod::track::eta) < cfgEtaCut) && (aod::track::pt > cfgMinPtCut) && (aod::track::pt < cfgMaxPtCut) && (aod::track::signed1Pt > cfgChargeCut);
  Partition<o2::aod::MyTracks> negative = (nabs(aod::track::eta) < cfgEtaCut) && (aod::track::pt > cfgMinPtCut) && (aod::track::pt < cfgMaxPtCut) && (aod::track::signed1Pt < cfgChargeCut);

  void init(o2::framework::InitContext&)
  {
    AxisSpec vtxZAxis = {100, -20, 20};
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    if (cfgUseCCDB) {
      ccdbApi.init(cfgCCDBURL); // Initializes ccdbApi when cfgUseCCDB is set to 'true'
    }
    pidModel = std::make_shared<PidONNXModel>(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, cfgTimestamp.value, cfgPid.value, static_cast<PidMLDetector>(cfgDetector.value), cfgCertainty.value);

    histos.add("hChargePos", ";z;", kTH1F, {{3, -1.5, 1.5}});
    histos.add("hChargeNeg", ";z;", kTH1F, {{3, -1.5, 1.5}});
    histos.add("hInvariantMass", ";M_{k^{+}k^{-}} (GeV/#it{c}^{2});", kTH1F, {{100, 0., 2.}});
    histos.add("hdEdXvsMomentum", ";P_{K^{+}K^{-}}; dE/dx in TPC (keV/cm)", kTH2F, {{100, 0., 4.}, {200, 20., 400.}});
  }

  void process(MyFilteredCollision const& coll, o2::aod::MyTracks const& /*tracks*/)
  {
    auto groupPositive = positive->sliceByCached(aod::track::collisionId, coll.globalIndex(), cache);
    auto groupNegative = negative->sliceByCached(aod::track::collisionId, coll.globalIndex(), cache);
    for (auto track : groupPositive) {
      histos.fill(HIST("hChargePos"), track.sign());
      if (pidModel.get()->applyModelBoolean(track)) {
        histos.fill(HIST("hdEdXvsMomentum"), track.p(), track.tpcSignal());
      }
    }

    for (auto track : groupNegative) {
      histos.fill(HIST("hChargeNeg"), track.sign());
      if (pidModel.get()->applyModelBoolean(track)) {
        histos.fill(HIST("hdEdXvsMomentum"), track.p(), track.tpcSignal());
      }
    }

    for (auto& [pos, neg] : combinations(soa::CombinationsFullIndexPolicy(groupPositive, groupNegative))) {
      if (!(pidModel.get()->applyModelBoolean(pos)) || !(pidModel.get()->applyModelBoolean(neg))) {
        continue;
      }

      TLorentzVector part1Vec;
      TLorentzVector part2Vec;
      float mMassOne = TDatabasePDG::Instance()->GetParticle(cfgPid.value)->Mass();
      float mMassTwo = TDatabasePDG::Instance()->GetParticle(cfgPid.value)->Mass();

      part1Vec.SetPtEtaPhiM(pos.pt(), pos.eta(), pos.phi(), mMassOne);
      part2Vec.SetPtEtaPhiM(neg.pt(), neg.eta(), neg.phi(), mMassTwo);

      TLorentzVector sumVec(part1Vec);
      sumVec += part2Vec;

      histos.fill(HIST("hInvariantMass"), sumVec.M());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<KaonPidTask>(cfgc)};
  return workflow;
}
