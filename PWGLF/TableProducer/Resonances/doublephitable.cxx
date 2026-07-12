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

/// \file doublephitable.cxx
/// \brief Selection of events with triplets and pairs for femtoscopic studies
///
/// \author Sourav Kundu, sourav.kundu@cern.ch

#include "PWGLF/DataModel/ReducedDoublePhiTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TMath.h>

#include <cstdint>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct doublephitable {

  // Produce derived tables
  Produces<aod::RedPhiEvents> redPhiEvents;
  Produces<aod::PhiTracks> phiTrack;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  /// Event selection
  struct : ConfigurableGroup {
    std::string prefix = std::string("EventSel");
    Configurable<bool> ConfEvtSelectZvtx{"ConfEvtSelectZvtx", true, "Event selection includes max. z-Vertex"};
    Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
    Configurable<bool> ConfEvtSel8{"ConfEvtSel8", true, "Event selection sel8"};
  } EventSel;

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 0.0f, "Accepted maximum Centrality"};
  // Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 100.0f, "Accepted minimum Centrality"};
  // track
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "ccdb-url"};
  Configurable<bool> useTrigger{"useTrigger", true, "use Trigger"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0.0, "cut TOF beta"};
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", -2.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isDeepAngle{"isDeepAngle", true, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {120, 0.98, 1.1}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  Configurable<float> minPhiMass{"minPhiMass", 1.01, "Minimum phi mass"};
  Configurable<float> maxPhiMass{"maxPhiMass", 1.03, "Maximum phi mass"};
  Configurable<float> nsigmaCutTPCPreSel{"nsigmacutTPCPreSel", 3.0, "Value of the TPC Nsigma cut Pre selection"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  // Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter PIDcutFilter = nabs(aod::pidtpc::tpcNSigmaKa) < nsigmaCutTPCPreSel;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  SliceCache cache;
  Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", ";; Number of events", 4, 0.0f, 4.0f)};
  // Histogram
  HistogramRegistry qaRegistry{"QAHistos", {
                                             {"hEventstat", "hEventstat", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
                                             {"hInvMassPhi", "hInvMassPhi", {HistType::kTH2F, {{40, 1.0f, 1.04f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtkaonTPC", "hNsigmaPtkaonTPC", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtkaonTOF", "hNsigmaPtkaonTOF", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                           },
                               OutputObjHandlingPolicy::AnalysisObject};

  double massKa = o2::constants::physics::MassKPlus;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(url.value);
    ccdbApi.init(url);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    zorroSummary.setObject(zorro.getZorroSummary());
    hProcessedEvents->GetXaxis()->SetBinLabel(1, "All Trigger events");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "Events with Double Phi Trigger");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "Events eith trigger and Evsel");
    hProcessedEvents->GetXaxis()->SetBinLabel(4, "Final Event");
  }
  template <typename T>
  bool isSelectedEvent(T const& col)
  {
    if (EventSel.ConfEvtSelectZvtx && std::abs(col.posZ()) > EventSel.ConfEvtZvtx) {
      return false;
    }
    if (EventSel.ConfEvtSel8 && !col.sel8()) {
      return false;
    }
    return true;
  }
  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (useGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsCrossedRows() > cfgTPCcluster)) {
      return false;
    }
    return true;
  }
  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (candidate.pt() < 0.5 && candidate.tpcNSigmaKa() > nsigmaCutTPC && candidate.tpcNSigmaKa() < 3.0) {
      return true;
    }
    if (candidate.pt() >= 0.5) {
      if (!candidate.hasTOF() && candidate.tpcNSigmaKa() > nsigmaCutTPC && candidate.tpcNSigmaKa() < 2.0) {
        return true;
      }
      if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && TMath::Sqrt(candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < nsigmaCutTOF) {
        return true;
      }
    }
    return false;
  }

  int currentRunNumber = -999;
  int lastRunNumber = -999;
  ROOT::Math::PxPyPzMVector KaonPlus, KaonMinus, PhiMesonMother, PhiVectorDummy, Phid1dummy, Phid2dummy;
  void processPhiReducedTable(EventCandidates::iterator const& collision, TrackCandidates const&, aod::BCsWithTimestamps const&)
  {
    o2::aod::ITSResponse itsResponse;
    bool keepEventDoublePhi = false;
    int numberPhi = 0;

    std::vector<int64_t> Phid1Index = {};
    std::vector<int64_t> Phid2Index = {};

    std::vector<float> Phid1Charge = {};
    std::vector<float> Phid2Charge = {};

    std::vector<float> Phid1TPC = {};
    std::vector<float> Phid2TPC = {};

    std::vector<float> Phid1TOF = {};
    std::vector<float> Phid2TOF = {};

    std::vector<int> Phid1TOFHit = {};
    std::vector<int> Phid2TOFHit = {};

    std::vector<ROOT::Math::PtEtaPhiMVector> phiresonance, phiresonanced1, phiresonanced2;

    int Npostrack = 0;
    int Nnegtrack = 0;
    float centrality = collision.centFT0M();
    currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    hProcessedEvents->Fill(0.5);
    bool zorroSelected = false;
    if (currentRunNumber != lastRunNumber) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fTriggerEventDoublePhi");
      zorro.populateHistRegistry(qaRegistry, bc.runNumber());
      lastRunNumber = currentRunNumber;
    }
    if (useTrigger) {
      zorroSelected = zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC());
    } else {
      zorroSelected = true;
    }
    if (zorroSelected) {
      hProcessedEvents->Fill(1.5);
    }
    if (zorroSelected && isSelectedEvent(collision)) {
      hProcessedEvents->Fill(2.5);
      auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      for (auto track1 : posThisColl) {
        // track selection
        if (!selectionTrack(track1)) {
          continue;
        }
        // PID check
        if (!selectionPID(track1)) {
          continue;
        }
        if (track1.pt() > 0.4 && track1.pt() < 1.0 && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1) > -2.0 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1) < 3.0)) {
          continue;
        }
        Npostrack = Npostrack + 1;
        qaRegistry.fill(HIST("hNsigmaPtkaonTPC"), track1.tpcNSigmaKa(), track1.pt());
        if (track1.hasTOF()) {
          qaRegistry.fill(HIST("hNsigmaPtkaonTOF"), track1.tofNSigmaKa(), track1.pt());
        }
        auto track1ID = track1.globalIndex();

        for (auto track2 : negThisColl) {
          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID) {
            continue;
          }
          // track selection
          if (!selectionTrack(track2)) {
            continue;
          }
          // PID check
          if (!selectionPID(track2)) {
            continue;
          }
          if (track2.pt() > 0.4 && track2.pt() < 1.0 && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) > -2.0 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) < 3.0)) {
            continue;
          }
          if (Npostrack == 1) {
            Nnegtrack = Nnegtrack + 1;
          }
          KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
          KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
          PhiMesonMother = KaonPlus + KaonMinus;
          if (PhiMesonMother.M() > minPhiMass && PhiMesonMother.M() < maxPhiMass) {
            numberPhi = numberPhi + 1;
            ROOT::Math::PtEtaPhiMVector temp1(track1.pt(), track1.eta(), track1.phi(), massKa);
            ROOT::Math::PtEtaPhiMVector temp2(track2.pt(), track2.eta(), track2.phi(), massKa);
            ROOT::Math::PtEtaPhiMVector temp3(PhiMesonMother.pt(), PhiMesonMother.eta(), PhiMesonMother.phi(), PhiMesonMother.M());

            phiresonanced1.push_back(temp1);
            phiresonanced2.push_back(temp2);
            phiresonance.push_back(temp3);

            Phid1Index.push_back(track1.globalIndex());
            Phid2Index.push_back(track2.globalIndex());

            Phid1Charge.push_back(track1.sign());
            Phid2Charge.push_back(track2.sign());

            Phid1TPC.push_back(track1.tpcNSigmaKa());
            Phid2TPC.push_back(track2.tpcNSigmaKa());
            auto d1TOFHit = -1;
            auto d2TOFHit = -1;
            auto d1TOF = -999.0;
            auto d2TOF = -999.0;
            if (track1.hasTOF()) {
              d1TOFHit = 1;
              d1TOF = track1.tofNSigmaKa();
            }
            if (track2.hasTOF()) {
              d2TOFHit = 1;
              d2TOF = track2.tofNSigmaKa();
            }
            Phid1TOF.push_back(d1TOF);
            Phid2TOF.push_back(d2TOF);
            Phid1TOFHit.push_back(d1TOFHit);
            Phid2TOFHit.push_back(d2TOFHit);
            qaRegistry.fill(HIST("hInvMassPhi"), PhiMesonMother.M(), PhiMesonMother.Pt());
          }
        }
      }
    } // select collision
    if (numberPhi > 1 && zorroSelected && Nnegtrack > 1 && Npostrack > 1) {
      keepEventDoublePhi = true;
      hProcessedEvents->Fill(3.5);
    }
    qaRegistry.fill(HIST("hEventstat"), 0.5);
    if (keepEventDoublePhi && numberPhi > 1 && (phiresonance.size() == phiresonanced1.size()) && (phiresonance.size() == phiresonanced2.size())) {
      qaRegistry.fill(HIST("hEventstat"), 1.5);
      /////////// Fill collision table///////////////
      redPhiEvents(bc.globalBC(), currentRunNumber, bc.timestamp(), collision.posZ(), collision.numContrib(), Npostrack, Nnegtrack, centrality);
      auto indexEvent = redPhiEvents.lastIndex();
      //// Fill track table for Phi//////////////////
      for (auto if1 = phiresonance.begin(); if1 != phiresonance.end(); ++if1) {
        auto i5 = std::distance(phiresonance.begin(), if1);
        PhiVectorDummy = phiresonance.at(i5);
        Phid1dummy = phiresonanced1.at(i5);
        Phid2dummy = phiresonanced2.at(i5);
        phiTrack(indexEvent, PhiVectorDummy.Px(), PhiVectorDummy.Py(), PhiVectorDummy.Pz(), Phid1dummy.Px(), Phid1dummy.Py(), Phid1dummy.Pz(), Phid2dummy.Px(), Phid2dummy.Py(), Phid2dummy.Pz(),
                 PhiVectorDummy.M(),
                 Phid1Index.at(i5), Phid2Index.at(i5),
                 Phid1Charge.at(i5), Phid2Charge.at(i5),
                 Phid1TPC.at(i5), Phid2TPC.at(i5), Phid1TOFHit.at(i5), Phid2TOFHit.at(i5), Phid1TOF.at(i5), Phid2TOF.at(i5));
      }
    }
  } // process
  PROCESS_SWITCH(doublephitable, processPhiReducedTable, "Process table creation for double phi", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<doublephitable>(cfg)};
}
