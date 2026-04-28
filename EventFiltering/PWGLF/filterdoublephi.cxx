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

/// \file filterdoublephi.cxx
/// \brief Selection of events with triplets and pairs for femtoscopic studies
///
/// \author Sourav Kundu, sourav.kundu@cern.ch

#include "../filterTables.h"

#include "PWGLF/DataModel/ReducedDoublePhiTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/BetheBlochAleph.h"
#include <Framework/Configurable.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>
#include <TDatabasePDG.h> //FIXME
#include <TMath.h>
#include <TPDGCode.h> // FIXME

#include <fairlogger/Logger.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct filterdoublephi {

  // Produce derived tables
  Produces<aod::DoublePhiFilters> tags;

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 0.0f, "Accepted maximum Centrality"};
  // Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 100.0f, "Accepted minimum Centrality"};
  // track
  Configurable<bool> isPtdepPID1{"isPtdepPID1", false, "use pt dep PID kplus"};
  Configurable<bool> isPtdepPID2{"isPtdepPID2", false, "use pt dep PID kminus"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0.0, "cut TOF beta"};
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", -2.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutTPCPreSel{"nsigmacutTPCPreSel", 3.0, "Value of the TPC Nsigma cut Pre selection"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isDeepAngle{"isDeepAngle", true, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {120, 0.98, 1.1}, "#it{M} (GeV/#it{c}^{2})"};
  ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
  Configurable<float> minPhiMass{"minPhiMass", 1.01, "Minimum phi mass"};
  Configurable<float> maxPhiMass{"maxPhiMass", 1.03, "Maximum phi mass"};
  Configurable<float> MinPhiPairPt{"MinPhiPairPt", 2.0, "Minimum phi pair Pt"};
  Configurable<float> MinPhiPairMass{"MinPhiPairMass", 2.5, "Minimum phi pair mass"};
  Configurable<float> MaxPhiPairMass{"MaxPhiPairMass", 3.0, "Max phi pair mass"};

  // Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  // Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter PIDcutFilter = nabs(aod::pidtpc::tpcNSigmaKa) < nsigmaCutTPCPreSel;

  // using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults>>;
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  SliceCache cache;
  Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  // Histogram
  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", ";; Number of events", 3, 0.0f, 3.0f)};
  HistogramRegistry qaRegistry{"QAHistos", {
                                             {"hInvMassPhi", "hInvMassPhi", {HistType::kTH2F, {{40, 1.0f, 1.04f}, {100, 0.0f, 10.0f}}}},
                                             {"hInvMassDoublePhi", "hInvMassDoublePhi", {HistType::kTH2F, {{1000, 2.0f, 3.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtkaonTPC", "hNsigmaPtkaonTPC", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtkaonTOF", "hNsigmaPtkaonTOF", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                           },
                               OutputObjHandlingPolicy::AnalysisObject};

  double massKa = o2::constants::physics::MassKPlus;

  void init(o2::framework::InitContext&)
  {
    hProcessedEvents->GetXaxis()->SetBinLabel(1, "All events");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "Events with double Phi without sel.");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, aod::filtering::TriggerEventDoublePhi::columnLabel());
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
  template <typename T>
  bool selectionPID2(const T& candidate)
  {
    if (candidate.pt() < 0.5 && candidate.tpcNSigmaKa() > nsigmaCutTPC && candidate.tpcNSigmaKa() < 3.0) {
      return true;
    }
    if (candidate.pt() >= 0.5 && candidate.pt() < 5.0) {
      if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && TMath::Sqrt(candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < nsigmaCutTOF) {
        return true;
      }
    }
    if (candidate.pt() >= 5.0 && candidate.tpcNSigmaKa() > nsigmaCutTPC && candidate.tpcNSigmaKa() < 2.0) {
      return true;
    }
    return false;
  }
  // deep angle cut on pair to remove photon conversion
  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.pt();
    pt2 = candidate2.pt();
    pz1 = candidate1.pz();
    pz2 = candidate2.pz();
    p1 = candidate1.p();
    p2 = candidate2.p();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (isDeepAngle && angle < cfgDeepAngle) {
      return false;
    }
    return true;
  }

  ROOT::Math::PxPyPzMVector KaonPlus, KaonMinus, PhiMesonMother, PhiVectorDummy, PhiVectorDummy2, PhiPair;
  void processPhiReducedTable(EventCandidates::iterator const& collision, TrackCandidates const&, aod::BCsWithTimestamps const&)
  {
    bool keepEventDoublePhi = false;
    int numberPhi = 0;
    o2::aod::ITSResponse itsResponse;
    std::vector<int64_t> Phid1Index = {};
    std::vector<int64_t> Phid2Index = {};
    std::vector<ROOT::Math::PtEtaPhiMVector> phiresonance, phiresonanced1, phiresonanced2;
    int Npostrack = 0;
    int Nnegtrack = 0;
    hProcessedEvents->Fill(0.5);
    if (collision.sel8()) {
      auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      for (auto track1 : posThisColl) {
        // track selection
        if (!selectionTrack(track1)) {
          continue;
        }
        // PID check
        if (isPtdepPID1 && !selectionPID2(track1)) {
          continue;
        }
        if (!isPtdepPID1 && !selectionPID(track1)) {
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
          // track selection
          if (!selectionTrack(track2)) {
            continue;
          }
          // PID check
          if (isPtdepPID2 && !selectionPID2(track2)) {
            continue;
          }
          if (!isPtdepPID2 && !selectionPID(track2)) {
            continue;
          }
          if (track2.pt() > 0.4 && track2.pt() < 1.0 && !(itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) > -2.0 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2) < 3.0)) {
            continue;
          }
          if (Npostrack == 1) {
            Nnegtrack = Nnegtrack + 1;
          }
          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID) {
            continue;
          }
          if (!selectionPair(track1, track2)) {
            continue;
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
            qaRegistry.fill(HIST("hInvMassPhi"), PhiMesonMother.M(), PhiMesonMother.Pt());
          }
        }
      }
    } // select collision
    if (numberPhi > 1 && Npostrack > 1 && Nnegtrack > 1 && (phiresonance.size() == phiresonanced1.size()) && (phiresonance.size() == phiresonanced2.size())) {
      hProcessedEvents->Fill(1.5);
      for (auto if1 = phiresonance.begin(); if1 != phiresonance.end(); ++if1) {
        auto i5 = std::distance(phiresonance.begin(), if1);
        PhiVectorDummy = phiresonance.at(i5);
        for (auto if2 = if1 + 1; if2 != phiresonance.end(); ++if2) {
          auto i6 = std::distance(phiresonance.begin(), if2);
          PhiVectorDummy2 = phiresonance.at(i6);
          PhiPair = PhiVectorDummy + PhiVectorDummy2;
          if (!(Phid1Index.at(i5) == Phid1Index.at(i6) || Phid2Index.at(i5) == Phid2Index.at(i6)) && PhiPair.M() > MinPhiPairMass && PhiPair.M() < MaxPhiPairMass && PhiPair.Pt() > MinPhiPairPt) {
            qaRegistry.fill(HIST("hInvMassDoublePhi"), PhiPair.M(), PhiPair.Pt());
            keepEventDoublePhi = true;
          }
        }
      }
    }
    if (keepEventDoublePhi) {
      hProcessedEvents->Fill(2.5);
    }
    tags(keepEventDoublePhi);
  } // process
  PROCESS_SWITCH(filterdoublephi, processPhiReducedTable, "Process table creation for double phi", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<filterdoublephi>(cfg, TaskName{"lf-doublephi-filter"})};
}
