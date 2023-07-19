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

/// \author Junlee Kim (jikim1290@gmail.com)

#include <TLorentzVector.h>

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "DataFormatsParameters/GRPObject.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct f0980analysis {
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> cfgMinPt{"cfgMinPt", 0.15, "Minimum transverse momentum for charged track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5, "Maximum transverse DCA"};
  Configurable<float> cfgMinDCAzToPVcut{"cfgMinDCAzToPVcut", 0.0, "Minimum longitudinal DCA"};
  Configurable<float> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0, "Maximum longitudinal DCA"};
  Configurable<float> cfgMaxTPCStandalone{"cfgMaxTPCStandalone", 2.0, "Maximum TPC PID as standalone"};
  Configurable<float> cfgMaxTPC{"cfgMaxTPC", 5.0, "Maximum TPC PID with TOF"};
  Configurable<float> cfgMaxTOF{"cfgMaxTOF", 3.0, "Maximum TOF PID with TPC"};
  Configurable<float> cfgMinRap{"cfgMinRap", -0.5, "Minimum rapidity for pair"};
  Configurable<float> cfgMaxRap{"cfgMaxRap", 0.5, "Maximum rapidity for pair"};
  Configurable<int> cfgMinTPCncr{"cfgMinTPCncr", 70, "minimum TPC cluster"};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> ptBinning = {
      0.0, 0.2, 0.4, 0.6, 0.8,
      1.0, 1.5, 2.0, 2.5, 3.0,
      3.5, 4.0, 4.5, 5.0, 6.0,
      7.0, 8.0, 10.0, 13.0, 20.0};

    AxisSpec centAxis = {20, 0, 100};
    AxisSpec ptAxis = {ptBinning};
    AxisSpec massAxis = {400, 0.2, 2.2};
    AxisSpec epAxis = {20, -constants::math::PI, constants::math::PI};

    AxisSpec PIDqaAxis = {120, -6, 6};
    AxisSpec pTqaAxis = {100, 0, 20};

    histos.add("hInvMass_f0980_US", "unlike invariant mass", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});
    histos.add("hInvMass_f0980_LSpp", "++ invariant mass", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});
    histos.add("hInvMass_f0980_LSmm", "-- invariant mass", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});

    histos.add("QA/Nsigma_TPC", "", {HistType::kTH2F, {pTqaAxis, PIDqaAxis}});
    histos.add("QA/Nsigma_TOF", "", {HistType::kTH2F, {pTqaAxis, PIDqaAxis}});
    histos.add("QA/TPC_TOF", "", {HistType::kTH2F, {PIDqaAxis, PIDqaAxis}});

    histos.print();
  }

  double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();

  template <typename TrackType>
  bool SelTrack(const TrackType track)
  {
    if (track.pt() < cfgMinPt)
      return false;
    if (std::fabs(track.eta()) > cfgMaxEta)
      return false;
    if (track.dcaXY() > cfgMaxDCArToPVcut)
      return false;
    if (track.dcaZ() < cfgMinDCAzToPVcut || track.dcaZ() > cfgMaxDCAzToPVcut)
      return false;

    return true;
  }

  template <typename TrackType>
  bool SelPion(const TrackType track)
  {
    if ((track.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) != aod::resodaughter::kHasTOF) {
      if (std::fabs(track.tpcNSigmaPi()) > cfgMaxTPCStandalone) {
        return false;
      }
    } else {
      if (std::fabs(track.tpcNSigmaPi()) > cfgMaxTPC || std::fabs(track.tofNSigmaPi()) > cfgMaxTOF) {
        return false;
      }
    }
    return true;
  }

  template <bool IsMC, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks)
  {
    TLorentzVector Pion1, Pion2, Reco;
    for (auto& [trk1, trk2] : combinations(CombinationsUpperIndexPolicy(dTracks, dTracks))) {

      if (trk1.index() == trk2.index()) {
        histos.fill(HIST("QA/Nsigma_TPC"), trk1.pt(), trk1.tpcNSigmaPi());
        histos.fill(HIST("QA/Nsigma_TOF"), trk1.pt(), trk1.tofNSigmaPi());
        histos.fill(HIST("QA/TPC_TOF"), trk1.tpcNSigmaPi(), trk1.tofNSigmaPi());
        continue;
      }

      if (!SelTrack(trk1) || !SelTrack(trk2))
        continue;
      if (!SelPion(trk1) || !SelPion(trk2))
        continue;

      Pion1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPi);
      Pion2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massPi);
      Reco = Pion1 + Pion2;

      if (Reco.Rapidity() > cfgMaxRap || Reco.Rapidity() < cfgMinRap)
        continue;

      if (trk1.sign() * trk2.sign() < 0) {
        histos.fill(HIST("hInvMass_f0980_US"), Reco.M(), Reco.Pt(), collision.multV0M());
      } else if (trk1.sign() > 0 && trk2.sign() > 0) {
        histos.fill(HIST("hInvMass_f0980_LSpp"), Reco.M(), Reco.Pt(), collision.multV0M());
      } else if (trk1.sign() < 0 && trk2.sign() < 0) {
        histos.fill(HIST("hInvMass_f0980_LSmm"), Reco.M(), Reco.Pt(), collision.multV0M());
      }
    }
  }

  void processData(aod::ResoCollision& collision,
                   aod::ResoTracks const& resotracks)
  {
    fillHistograms<false>(collision, resotracks);
  };
  PROCESS_SWITCH(f0980analysis, processData, "Process Event for data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<f0980analysis>(cfgc, TaskName{"lf-f0980analysis"})};
}
