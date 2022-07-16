#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/TrackSelectorPID.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include <CCDB/BasicCCDBManager.h>
#include "DataFormatsParameters/GRPObject.h"
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct k892sample {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  framework::Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB


  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    AxisSpec vtxZAxis = {100, -20, 20};

    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 10., 20.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    histos.add("EventQA/Centrality", "Centrality distribution (V0M)", kTH1F, {centAxis});
    histos.add("EventQA/VtxZBeforeSel", "Vertex distribution in Z;Z (cm)", kTH1F, {vtxZAxis});
    histos.add("EventQA/VtxZAfterSel", "Vertex distribution in Z;Z (cm)", kTH1F, {vtxZAxis});
    // Mass QA (quick check)
    histos.add("k892invmass", "Invariant mass of K(892)", kTH1F, {{500, 0.5, 1.1, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("antik892invmass", "Invariant mass of Anti-K(892)", kTH1F, {{500, 0.5, 1.1, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("k892invmassME", "Invariant mass of K(892) mixed event", kTH1F, {{500, 0.5, 1.1, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("antik892invmassME", "Invariant mass of Anti-K(892) mixed event", kTH1F, {{500, 0.5, 1.1, "Invariant Mass (GeV/#it{c}^2)"}});

    // 3d histogram
    histos.add("h3k892invmass", "Invariant mass of K(892)", kTH3F, {{100, 0.0f, 100.0f}, {100, 0.0f, 10.0f}, {500, 0.5, 1.1}});
    histos.add("h3antik892invmass", "Invariant mass of Anti-K(892)", kTH3F, {{100, 0.0f, 100.0f}, {100, 0.0f, 10.0f}, {500, 0.5, 1.1}});
  }

  double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
  double massKa = TDatabasePDG::Instance()->GetParticle(kKPlus)->Mass();

  void process(o2::aod::ResoCollision& inputCollision,
               aod::Reso2TrackTrackDatas const& reso2trktrk, aod::Reso2TracksPIDExt const& track)
  {

    histos.fill(HIST("EventQA/VtxZAfterSel"), inputCollision.posZ());

    // fill centrality histos
    histos.fill(HIST("EventQA/Centrality"), inputCollision.multV0M());

    for (auto reso2 : reso2trktrk) {
      auto trackPos = (reso2.trk1Sign() > 0) ? reso2.track1_as<aod::Reso2TracksPIDExt>() : reso2.track2_as<aod::Reso2TracksPIDExt>(); // positive daughter
      auto trackNeg = (reso2.trk1Sign() > 0) ? reso2.track2_as<aod::Reso2TracksPIDExt>() : reso2.track1_as<aod::Reso2TracksPIDExt>(); // negative daughter
      bool isK892 = false;
      bool isAntiK892 = false;

      // PID cuts
      if ((std::abs(trackPos.tpcNSigmaKa()) < 3) && (std::abs(trackNeg.tpcNSigmaPi()) < 3)) // pi- + K+
        isK892 = true;
      if ((std::abs(trackNeg.tpcNSigmaKa()) < 3) && (std::abs(trackPos.tpcNSigmaPi()) < 3)) // K- + pi+ (anti)
        isAntiK892 = true;
      if (!isK892 && !isAntiK892)
        continue;

      // TOF PID cut (if available)
      if (isK892 && trackPos.hasTOF() && (std::abs(trackPos.tofNSigmaKa()) > 3))
        continue;
      if (isAntiK892 && trackNeg.hasTOF() && (std::abs(trackNeg.tofNSigmaKa()) > 3))
        continue;

      auto arrMom = array{
        array{reso2.pxTrk1(), reso2.pyTrk1(), reso2.pzTrk1()},
        array{reso2.pxTrk2(), reso2.pyTrk2(), reso2.pzTrk2()}};
      auto arrMass = (isK892) ? array{massKa, massPi} : array{massPi, massKa};
      auto mass = RecoDecay::m(arrMom, arrMass);
      if (isK892) {
        histos.fill(HIST("k892invmass"), mass);
        histos.fill(HIST("h3k892invmass"), inputCollision.multV0M(), reso2.pt(), mass);
      }
      if (isAntiK892) { // Anti-matter
        histos.fill(HIST("antik892invmass"), mass);
        histos.fill(HIST("h3antik892invmass"), inputCollision.multV0M(), reso2.pt(), mass);
      }
    }
  }
  // Processing Event Mixing
  void processME(o2::aod::ResoCollisions& inputCollision,
                 o2::aod::BCsWithTimestamps const&, aod::Reso2TrackTrackDatas & reso2trktrk, aod::Reso2TracksPIDExt const& track)
  {
    BinningPolicy<aod::collision::PosZ, aod::resocollision::MultV0M> colBinning{{CfgVtxBins, CfgMultBins}, true};

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, inputCollision, inputCollision)) {

      auto magFieldTesla1 = 0.0;
      static o2::parameters::GRPObject* grpo1 = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", collision1.timestamp());
      if (!grpo1) {
        magFieldTesla1 = 0;
      } else {
        // LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", collision1.timestamp(), grpo1->getNominalL3Field());
        magFieldTesla1 = 0.1 * (grpo1->getNominalL3Field());
      }

      auto magFieldTesla2 = 0.0;
      static o2::parameters::GRPObject* grpo2 = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", collision2.timestamp());
      if (!grpo2) {
        magFieldTesla2 = 0;
      } else {
        // LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", collision2.timestamp(), grpo2->getNominalL3Field());
        magFieldTesla2 = 0.1 * (grpo1->getNominalL3Field());
      }

      auto reso2trktrkPartOne = reso2trktrk.sliceByCached(aod::reso2trktrkdata::collisionId, collision1.globalIndex());
      auto reso2trktrkPartTwo = reso2trktrk.sliceByCached(aod::reso2trktrkdata::collisionId, collision2.globalIndex());

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      for (auto& [reso1, reso2] : combinations(CombinationsStrictlyUpperIndexPolicy(reso2trktrkPartOne, reso2trktrkPartTwo))) {
        auto trackPos = (reso1.trk1Sign() > 0) ? reso1.track1_as<aod::Reso2TracksPIDExt>() : reso1.track2_as<aod::Reso2TracksPIDExt>(); // positive daughter
        auto trackNeg = (reso2.trk1Sign() > 0) ? reso2.track2_as<aod::Reso2TracksPIDExt>() : reso2.track1_as<aod::Reso2TracksPIDExt>(); // negative daughter
        bool isK892 = false;
        bool isAntiK892 = false;

        // PID cuts
        if ((std::abs(trackPos.tpcNSigmaKa()) < 3) && (std::abs(trackNeg.tpcNSigmaPi()) < 3)) // pi- + K+
          isK892 = true;
        if ((std::abs(trackNeg.tpcNSigmaKa()) < 3) && (std::abs(trackPos.tpcNSigmaPi()) < 3)) // K- + pi+ (anti)
          isAntiK892 = true;
        if (!isK892 && !isAntiK892)
          continue;

        // TOF PID cut (if available)
        if (isK892 && trackPos.hasTOF() && (std::abs(trackPos.tofNSigmaKa()) > 3))
          continue;
        if (isAntiK892 && trackNeg.hasTOF() && (std::abs(trackNeg.tofNSigmaKa()) > 3))
          continue;

        auto arrMom = array{
          array{reso1.pxTrk1(), reso1.pyTrk1(), reso1.pzTrk1()},
          array{reso2.pxTrk2(), reso2.pyTrk2(), reso2.pzTrk2()}};
        auto arrMass = (isK892) ? array{massKa, massPi} : array{massPi, massKa};
        auto mass = RecoDecay::m(arrMom, arrMass);
        if (isK892) {
          histos.fill(HIST("k892invmassME"), mass);
        }
        if (isAntiK892) { // Anti-matter
          histos.fill(HIST("antik892invmassME"), mass);
        }
      }
    }
  };
  PROCESS_SWITCH(k892sample, processME, "Process EventMixing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<k892sample>(cfgc)};
  return workflow;
}