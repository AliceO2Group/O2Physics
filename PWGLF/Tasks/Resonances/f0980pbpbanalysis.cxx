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

#include <cmath>
#include <array>
#include <cstdlib>
#include <chrono>
#include <string>

#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TVector2.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include <TMath.h>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"

#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"

#include "CommonConstants/PhysicsConstants.h"

#include "ReconstructionDataFormats/Track.h"

#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct f0980pbpbanalysis {
  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0, "PV selection"};
  Configurable<bool> cfgQvecSel{"cfgQvecSel", true, "Reject events when no QVector"};
  Configurable<bool> cfgOccupancySel{"cfgOccupancySe", false, "Occupancy selection"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", -100, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<bool> cfgNCollinTR{"cfgNCollinTR", false, "Additional selection for the number of coll in time range"};
  Configurable<bool> cfgPVSel{"cfgPVSel", false, "Additional PV selection flag for syst"};
  Configurable<float> cfgPV{"cfgPV", 8.0, "Additional PV selection range for syst"};

  Configurable<float> cfgCentSel{"cfgCentSel", 80., "Centrality selection"};
  Configurable<int> cfgCentEst{"cfgCentEst", 1, "Centrality estimator, 1: FT0C, 2: FT0M"};

  Configurable<float> cfgMinPt{"cfgMinPt", 0.15, "Minimum transverse momentum for charged track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5, "Maximum transverse DCA"};
  Configurable<float> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0, "Maximum longitudinal DCA"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<float> cfgRatioTPCRowsOverFindableCls{"cfgRatioTPCRowsOverFindableCls", 0.8, "TPC Crossed Rows to Findable Clusters"};

  Configurable<float> cfgMinRap{"cfgMinRap", -0.5, "Minimum rapidity for pair"};
  Configurable<float> cfgMaxRap{"cfgMaxRap", 0.5, "Maximum rapidity for pair"};

  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};

  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"}; // TOF
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 5.0, "TPC nSigma cut for Pion"}; // TPC
  Configurable<double> cMaxTPCnSigmaPionS{"cMaxTPCnSigmaPionS", 3.0, "TPC nSigma cut for Pion as a standalone"};
  Configurable<bool> cfgUSETOF{"cfgUSETOF", false, "TPC usage"};

  Configurable<int> cfgnMods{"cfgnMods", 1, "The number of modulations of interest starting from 2"};
  Configurable<int> cfgNQvec{"cfgNQvec", 7, "The number of total Qvectors for looping over the task"};

  Configurable<std::string> cfgQvecDetName{"cfgQvecDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgQvecRefAName{"cfgQvecRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgQvecRefBName{"cfgQvecRefBName", "TPCneg", "The name of detector for reference B"};

  ConfigurableAxis massAxis{"massAxis", {400, 0.2, 2.2}, "Invariant mass axis"};
  ConfigurableAxis ptAxis{"ptAxis", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 13.0, 20.0}, "Transverse momentum Binning"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100}, "Centrality interval"};

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;

  int DetId;
  int RefAId;
  int RefBId;

  int QvecDetInd;
  int QvecRefAInd;
  int QvecRefBInd;

  float centrality;

  double angle;
  double relPhi;

  double massPi = o2::constants::physics::MassPionCharged;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgMaxEta && nabs(aod::track::pt) > cfgMinPt);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgMaxDCArToPVcut) && (nabs(aod::track::dcaZ) < cfgMaxDCAzToPVcut);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults, aod::Qvectors>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi>>;

  template <typename T>
  int GetDetId(const T& name)
  {
    if (name.value == "FT0C") {
      return 0;
    } else if (name.value == "FT0A") {
      return 1;
    } else if (name.value == "FT0M") {
      return 2;
    } else if (name.value == "FV0A") {
      return 3;
    } else if (name.value == "TPCpos") {
      return 4;
    } else if (name.value == "TPCneg") {
      return 5;
    } else {
      return 0;
    }
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision)
  {
    if (!collision.sel8()) {
      return 0;
    }

    if (cfgCentSel < centrality) {
      return 0;
    }
    /*
        auto multNTracksPV = collision.multNTracksPV();
        if (multNTracksPV < fMultPVCutLow->Eval(centrality)) {
          return 0;
        }
        if (multNTracksPV > fMultPVCutHigh->Eval(centrality)) {
          return 0;
        }
    */
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    if (cfgQvecSel && (collision.qvecAmp()[DetId] < 1e-4 || collision.qvecAmp()[RefAId] < 1e-4 || collision.qvecAmp()[RefAId] < 1e-4)) {
      return 0;
    }
    if (cfgOccupancySel && (collision.trackOccupancyInTimeRange() > cfgMaxOccupancy || collision.trackOccupancyInTimeRange() < cfgMinOccupancy)) {
      return 0;
    }
    if (cfgNCollinTR && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
    }
    if (cfgPVSel && std::abs(collision.posZ()) > cfgPV) {
      return 0;
    }

    return 1;
  } // event selection

  template <typename TrackType>
  bool trackSelected(const TrackType track)
  {
    if (std::abs(track.pt()) < cfgMinPt) {
      return 0;
    }
    if (std::fabs(track.eta()) > cfgMaxEta) {
      return 0;
    }
    if (std::fabs(track.dcaXY()) > cfgMaxDCArToPVcut) {
      return 0;
    }
    if (std::fabs(track.dcaZ()) > cfgMaxDCAzToPVcut) {
      return 0;
    }
    if (track.tpcNClsFound() < cfgTPCcluster) {
      return 0;
    }
    if (cfgPVContributor && !track.isPVContributor()) {
      return 0;
    }
    if (cfgPrimaryTrack && !track.isPrimaryTrack()) {
      return 0;
    }
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA()) {
      return 0;
    }
    if (track.tpcCrossedRowsOverFindableCls() < cfgRatioTPCRowsOverFindableCls) {
      return 0;
    }

    return 1;
  }

  template <typename TrackType>
  bool PIDSelected(const TrackType track)
  {
    if (cfgUSETOF) {
      if (std::fabs(track.tofNSigmaPi()) > cMaxTOFnSigmaPion) {
        return 0;
      }
      if (std::fabs(track.tpcNSigmaPi()) > cMaxTPCnSigmaPion) {
        return 0;
      }
    }
    if (std::fabs(track.tpcNSigmaPi()) > cMaxTPCnSigmaPionS) {
      return 0;
    }

    return 1;
  }

  template <bool IsMC, typename CollisionType, typename TracksType>
  void FillHistograms(const CollisionType& collision,
                      const TracksType& dTracks, int nmode)
  {
    QvecDetInd = DetId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefAInd = RefAId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    QvecRefBInd = RefBId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;

    double eventPlaneDet = TMath::ATan2(collision.qvecIm()[QvecDetInd], collision.qvecRe()[QvecDetInd]) / static_cast<float>(nmode);
    double eventPlaneRefA = TMath::ATan2(collision.qvecIm()[QvecRefAInd], collision.qvecRe()[QvecRefAInd]) / static_cast<float>(nmode);
    double eventPlaneRefB = TMath::ATan2(collision.qvecIm()[QvecRefBInd], collision.qvecRe()[QvecRefBInd]) / static_cast<float>(nmode);

    histos.fill(HIST("QA/EPhist"), centrality, eventPlaneDet);
    histos.fill(HIST("QA/EPResAB"), centrality, TMath::Cos(static_cast<float>(nmode) * (eventPlaneDet - eventPlaneRefA)));
    histos.fill(HIST("QA/EPResAC"), centrality, TMath::Cos(static_cast<float>(nmode) * (eventPlaneDet - eventPlaneRefB)));
    histos.fill(HIST("QA/EPResBC"), centrality, TMath::Cos(static_cast<float>(nmode) * (eventPlaneRefA - eventPlaneRefB)));

    TLorentzVector Pion1, Pion2, Reco;
    for (auto& [trk1, trk2] :
         combinations(CombinationsUpperIndexPolicy(dTracks, dTracks))) {
      if (trk1.index() == trk2.index()) {
        if (!trackSelected(trk1))
          continue;
        histos.fill(HIST("QA/Nsigma_TPC"), trk1.pt(), trk1.tpcNSigmaPi());
        histos.fill(HIST("QA/Nsigma_TOF"), trk1.pt(), trk1.tofNSigmaPi());
        histos.fill(HIST("QA/TPC_TOF"), trk1.tpcNSigmaPi(), trk1.tofNSigmaPi());
        continue;
      }

      if (!trackSelected(trk1) || !trackSelected(trk2))
        continue;
      if (!PIDSelected(trk1) || !PIDSelected(trk2))
        continue;

      if (trk1.index() == trk2.index()) {
        histos.fill(HIST("QA/Nsigma_TPC_selected"), trk1.pt(), trk1.tpcNSigmaPi());
        histos.fill(HIST("QA/Nsigma_TOF_selected"), trk1.pt(), trk1.tofNSigmaPi());
        histos.fill(HIST("QA/TPC_TOF_selected"), trk1.tpcNSigmaPi(), trk1.tofNSigmaPi());
      }

      Pion1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPi);
      Pion2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massPi);
      Reco = Pion1 + Pion2;

      if (Reco.Rapidity() > cfgMaxRap || Reco.Rapidity() < cfgMinRap)
        continue;

      relPhi = TVector2::Phi_0_2pi((Reco.Phi() - eventPlaneDet) * static_cast<float>(nmode));

      if (trk1.sign() * trk2.sign() < 0) {
        histos.fill(HIST("hInvMass_f0980_US_EPA"), Reco.M(), Reco.Pt(), centrality, relPhi);
        /*
                if constexpr (IsMC) {
                  if (abs(trk1.pdgCode()) != 211 || abs(trk2.pdgCode()) != 211)
                    continue;
                  if (trk1.motherId() != trk2.motherId())
                    continue;
                  if (abs(trk1.motherPDG()) != 9010221)
                    continue;
                  histos.fill(HIST("MCL/hpT_f0980_REC"), Reco.M(), Reco.Pt(), centrality);
                }
        */
      } else if (trk1.sign() > 0 && trk2.sign() > 0) {
        histos.fill(HIST("hInvMass_f0980_LSpp_EPA"), Reco.M(), Reco.Pt(), centrality, relPhi);
      } else if (trk1.sign() < 0 && trk2.sign() < 0) {
        histos.fill(HIST("hInvMass_f0980_LSmm_EPA"), Reco.M(), Reco.Pt(), centrality, relPhi);
      }
    }
  }

  void init(o2::framework::InitContext&)
  {
    AxisSpec epAxis = {6, 0.0, 2.0 * constants::math::PI};
    AxisSpec centQaAxis = {110, 0, 110};
    AxisSpec vzQaAxis = {100, -20, 20};
    AxisSpec PIDqaAxis = {100, -10, 10};
    AxisSpec pTqaAxis = {200, 0, 20};
    AxisSpec epQaAxis = {100, -1.0 * constants::math::PI, constants::math::PI};
    AxisSpec epresAxis = {102, -1.02, 1.02};

    histos.add("QA/CentDist", "", {HistType::kTH1F, {centQaAxis}});
    histos.add("QA/Vz", "", {HistType::kTH1F, {vzQaAxis}});

    histos.add("QA/Nsigma_TPC", "", {HistType::kTH2F, {pTqaAxis, PIDqaAxis}});
    histos.add("QA/Nsigma_TOF", "", {HistType::kTH2F, {pTqaAxis, PIDqaAxis}});
    histos.add("QA/TPC_TOF", "", {HistType::kTH2F, {PIDqaAxis, PIDqaAxis}});

    histos.add("QA/Nsigma_TPC_selected", "", {HistType::kTH2F, {pTqaAxis, PIDqaAxis}});
    histos.add("QA/Nsigma_TOF_selected", "", {HistType::kTH2F, {pTqaAxis, PIDqaAxis}});
    histos.add("QA/TPC_TOF_selected", "", {HistType::kTH2F, {PIDqaAxis, PIDqaAxis}});

    histos.add("QA/EPhist", "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
    histos.add("QA/EPResAB", "", {HistType::kTH2F, {centQaAxis, epresAxis}});
    histos.add("QA/EPResAC", "", {HistType::kTH2F, {centQaAxis, epresAxis}});
    histos.add("QA/EPResBC", "", {HistType::kTH2F, {centQaAxis, epresAxis}});

    histos.add("hInvMass_f0980_US_EPA", "unlike invariant mass",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, epAxis}});
    histos.add("hInvMass_f0980_LSpp_EPA", "++ invariant mass",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, epAxis}});
    histos.add("hInvMass_f0980_LSmm_EPA", "-- invariant mass",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, epAxis}});

    //    if (doprocessMCLight) {
    //      histos.add("MCL/hpT_f0980_GEN", "generated f0 signals", HistType::kTH1F, {pTqaAxis});
    //      histos.add("MCL/hpT_f0980_REC", "reconstructed f0 signals", HistType::kTH3F, {massAxis, pTqaAxis, centAxis});
    //    }

    DetId = GetDetId(cfgQvecDetName);
    RefAId = GetDetId(cfgQvecRefAName);
    RefBId = GetDetId(cfgQvecRefBName);

    if (DetId == RefAId || DetId == RefBId || RefAId == RefBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      DetId = 0;
      RefAId = 4;
      RefBId = 5;
    }

    fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
    fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);

    ccdb->setURL(cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  void processData(EventCandidates::iterator const& collision,
                   TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    if (cfgCentEst == 1) {
      centrality = collision.centFT0C();
    } else if (cfgCentEst == 2) {
      centrality = collision.centFT0M();
    }
    if (!eventSelected(collision)) {
      return;
    }
    histos.fill(HIST("QA/CentDist"), centrality, 1.0);
    histos.fill(HIST("QA/Vz"), collision.posZ(), 1.0);

    FillHistograms<false>(collision, tracks, 2); // second order
  };
  PROCESS_SWITCH(f0980pbpbanalysis, processData, "Process Event for data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<f0980pbpbanalysis>(cfgc, TaskName{"lf-f0980pbpbanalysis"})};
}
