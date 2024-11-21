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

// C++/ROOT includes.
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TMath.h>
#include <TVector2.h>

// o2Physics includes.
#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"

#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "CommonConstants/PhysicsConstants.h"

// o2 includes.

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace CMETrackPIDcolums
{
DECLARE_SOA_COLUMN(NPIDFlag, Npid, int8_t); // Flag tracks without proper binning as -1, and indicate type of particle 0->un-Id, 1->pion, 2->kaon, 3->proton
} // namespace CMETrackPIDcolums
DECLARE_SOA_TABLE(Flags, "AOD", "Flags", CMETrackPIDcolums::NPIDFlag);
} // namespace o2::aod

using TracksPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>;
struct FillPIDcolums {
  Configurable<float> cfgnSigmaCutTPCPi{"cfgnSigmaCutTPCPi", 3.0, "Value of the TPC Nsigma cut for pions"};
  Configurable<float> cfgnSigmaCutTPCKa{"cfgnSigmaCutTPCKa", 3.0, "Value of the TPC Nsigma cut for kaons"};
  Configurable<float> cfgnSigmaCutTPCPr{"cfgnSigmaCutTPCPr", 3.0, "Value of the TPC Nsigma cut for protons"};
  Configurable<float> cfgnSigmaCutTOFPi{"cfgnSigmaCutTOFPi", 3.0, "Value of the TOF Nsigma cut for pions"};
  Configurable<float> cfgnSigmaCutTOFKa{"cfgnSigmaCutTOFKa", 3.0, "Value of the TOF Nsigma cut for kaons"};
  Configurable<float> cfgnSigmaCutTOFPr{"cfgnSigmaCutTOFPr", 3.0, "Value of the TOF Nsigma cut for protons"};
  Configurable<float> cfgnSigmaCutCombine{"cfgnSigmaCutCombine", 3.0, "Value of the Combined Nsigma cut"};
  Configurable<float> cfgPtMaxforTPCOnlyPID{"cfgPtMaxforTPCOnlyPID", 0.4, "Maxmium track pt for TPC only PID,only when onlyTOF and onlyTOFHIT closed"};
  Configurable<float> cfgMinPtPID{"cfgMinPtPID", 0.15, "Minimum track #P_{t} for PID"};
  Configurable<float> cfgMaxEtaPID{"cfgMaxEtaPID", 0.8, "Maximum track #eta for PID"};

  ConfigurableAxis cfgrigidityBins{"cfgrigidityBins", {200, -10.f, 10.f}, "Binning for rigidity #it{p}^{TPC}/#it{z}"};
  ConfigurableAxis cfgdedxBins{"cfgdedxBins", {1000, 0.f, 1000.f}, "Binning for dE/dx"};
  ConfigurableAxis cfgnSigmaBins{"cfgnSigmaBins", {200, -5.f, 5.f}, "Binning for n sigma"};
  ConfigurableAxis cfgaxisptPID{"cfgaxisptPID", {24, 0, 12}, ""};

  Configurable<bool> onlyTOF{"onlyTOF", false, "only TOF tracks"};
  Configurable<bool> onlyTOFHIT{"onlyTOFHIT", false, "accept only TOF hit tracks at high pt"};
  bool onlyTPC = true;

  template <typename TrackType>
  bool SelTrack_PID(const TrackType track)
  {
    if (!(track.pt() > cfgMinPtPID))
      return false;
    if (!(std::abs(track.eta()) < cfgMaxEtaPID))
      return false;
    if (!track.passedITSNCls())
      return false;
    if (!track.passedITSChi2NDF())
      return false;
    if (!track.passedITSHits())
      return false;
    if (!track.passedTPCCrossedRowsOverNCls())
      return false;
    if (!track.passedTPCChi2NDF())
      return false;
    if (!track.passedDCAxy())
      return false;
    if (!track.passedDCAz())
      return false;
    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate, int8_t PID)
  {
    if (candidate.pt() > cfgPtMaxforTPCOnlyPID) {
      onlyTPC = false;
    }

    if (PID == 0) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < cfgnSigmaCutTOFPi) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < cfgnSigmaCutTOFPi) {
          return true;
        }
        if (!candidate.hasTOF() &&
            std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPCPi) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPCPi) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (cfgnSigmaCutCombine * cfgnSigmaCutCombine)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPCPi) {
          return true;
        }
      }
    } else if (PID == 1) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < cfgnSigmaCutTOFKa) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < cfgnSigmaCutTOFKa) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPCPi) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPCPi) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (cfgnSigmaCutCombine * cfgnSigmaCutCombine)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPCPi) {
          return true;
        }
      }
    } else if (PID == 2) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < cfgnSigmaCutTOFPr) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < cfgnSigmaCutTOFPr) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPCPr) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPCPr) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < (cfgnSigmaCutCombine * cfgnSigmaCutCombine)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPCPr) {
          return true;
        }
      }
    }
    return false;
  }

  HistogramRegistry histosQA{"histosQAPID", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    AxisSpec axisRigidity{cfgrigidityBins, "#it{p}^{TPC}/#it{z}"};
    AxisSpec axisdEdx{cfgdedxBins, "d#it{E}/d#it{x}"};
    AxisSpec axisnSigma{cfgnSigmaBins, "n_{#sigma}TPC"};
    AxisSpec axisPtPID{cfgaxisptPID, "#it{p}_{T}"};

    histosQA.add(Form("QA/PID/histdEdxTPC_All"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
    histosQA.add(Form("QA/PID/histdEdxTPC_Pi"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
    histosQA.add(Form("QA/PID/histnSigma_Pi"), "", {HistType::kTH1F, {axisnSigma}});
    histosQA.add(Form("QA/PID/histnSigma_Pt_Pi"), "", {HistType::kTH2F, {axisPtPID, axisnSigma}});
    histosQA.add(Form("QA/PID/histdEdxTPC_Ka"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
    histosQA.add(Form("QA/PID/histnSigma_Ka"), "", {HistType::kTH1F, {axisnSigma}});
    histosQA.add(Form("QA/PID/histnSigma_Pt_Ka"), "", {HistType::kTH2F, {axisPtPID, axisnSigma}});
    histosQA.add(Form("QA/PID/histdEdxTPC_Pr"), "", {HistType::kTH2F, {axisRigidity, axisdEdx}});
    histosQA.add(Form("QA/PID/histnSigma_Pr"), "", {HistType::kTH1F, {axisnSigma}});
    histosQA.add(Form("QA/PID/histnSigma_Pt_Pr"), "", {HistType::kTH2F, {axisPtPID, axisnSigma}});
  }
  Produces<aod::Flags> PIDCMEtable;
  void process(TracksPID const& tracks)
  {
    int8_t PID_flag;
    for (auto& track : tracks) {
      if (!SelTrack_PID(track)) {
        PID_flag = -1;
      } else {
        histosQA.fill(HIST("QA/PID/histdEdxTPC_All"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
        float nsigma_array[3] = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
        PID_flag = 0;
        for (int8_t i = 0; i < 3; i++) {
          if (selectionPID(track, i))
            PID_flag = PID_flag * 10 + i + 1;
          if (PID_flag > 10) {                                                                     // If a track is identified as two different tracks.
            if (fabs(nsigma_array[(PID_flag / 10) - 1]) < fabs(nsigma_array[(PID_flag % 10) - 1])) // The track is identified as the particle whose |nsigma| is the least.
              PID_flag /= 10;
            else
              PID_flag %= 10;
          }
        }

        switch (PID_flag) {
          case 1:
            histosQA.fill(HIST("QA/PID/histdEdxTPC_Pi"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
            histosQA.fill(HIST("QA/PID/histnSigma_Pi"), track.tpcNSigmaPi());
            histosQA.fill(HIST("QA/PID/histnSigma_Pt_Pi"), track.pt(), track.tpcNSigmaPi());
            break;
          case 2:
            histosQA.fill(HIST("QA/PID/histdEdxTPC_Ka"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
            histosQA.fill(HIST("QA/PID/histnSigma_Ka"), track.tpcNSigmaKa());
            histosQA.fill(HIST("QA/PID/histnSigma_Pt_Ka"), track.pt(), track.tpcNSigmaKa());
            break;
          case 3:
            histosQA.fill(HIST("QA/PID/histdEdxTPC_Pr"), track.sign() * track.tpcInnerParam(), track.tpcSignal());
            histosQA.fill(HIST("QA/PID/histnSigma_Pr"), track.tpcNSigmaPr());
            histosQA.fill(HIST("QA/PID/histnSigma_Pt_Pr"), track.pt(), track.tpcNSigmaPr());
            break;
        }
      }
      PIDCMEtable(PID_flag);
    }
  }
};

struct pidcme {
  HistogramRegistry histosQA{"histosmain", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::vector<int>> cfgnMods{"cfgnMods", {2}, "Modulation of interest"};
  Configurable<std::string> cfgDetName{"cfgDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "TPCneg", "The name of detector for reference B"};

  Configurable<int> cfgnTotalSystem{"cfgnTotalSystem", 7, "total qvector number"};

  Configurable<float> cfgMinPt{"cfgMinPt", 0.15, "Minimum transverse momentum for charged track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.1, "Maximum transverse DCA"};
  Configurable<float> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 1.0, "Maximum longitudinal DCA"};

  ConfigurableAxis cfgaxisQvecF{"cfgaxisQvecF", {300, -1, 1}, ""};
  ConfigurableAxis cfgaxisQvec{"cfgaxisQvec", {100, -3, 3}, ""};
  ConfigurableAxis cfgaxisCent{"cfgaxisCent", {100, 0, 100}, ""};

  ConfigurableAxis cfgaxiscos{"cfgaxiscos", {102, -1.02, 1.02}, ""};
  ConfigurableAxis cfgaxispt{"cfgaxispt", {100, 0, 10}, ""};
  ConfigurableAxis cfgaxisCentMerged{"cfgaxisCentMerged", {20, 0, 100}, ""};

  ConfigurableAxis cfgaxissumpt{"cfgaxissumpt", {7, 1, 8}, "Binning for #gamma and #delta pt(particle1 + particle2)"};
  ConfigurableAxis cfgaxisdeltaeta{"cfgaxisdeltaeta", {5, 0, 1}, "Binning for #gamma and #delta |#eta(particle1 - particle2)|"};

  Configurable<bool> OpenCME = {"cfgkOpeanCME", true, "open PID CME"};

  EventPlaneHelper helperEP;
  SliceCache cache;

  unsigned int mult1, mult2, mult3;
  int DetId;
  int RefAId;
  int RefBId;

  template <typename T>
  int GetDetId(const T& name)
  {
    if (name.value == "BPos" || name.value == "BNeg" || name.value == "BTot") {
      LOGF(warning, "Using deprecated label: %s. Please use TPCpos, TPCneg, TPCall instead.", name.value);
    }
    if (name.value == "FT0C") {
      return 0;
    } else if (name.value == "FT0A") {
      return 1;
    } else if (name.value == "FT0M") {
      return 2;
    } else if (name.value == "FV0A") {
      return 3;
    } else if (name.value == "TPCpos" || name.value == "BPos") {
      return 4;
    } else if (name.value == "TPCneg" || name.value == "BNeg") {
      return 5;
    } else if (name.value == "TPCall" || name.value == "BTot") {
      return 6;
    } else {
      return 0;
    }
  }

  Filter col = aod::evsel::sel8 == true;
  Filter collisionFilter = (nabs(aod::collision::posZ) < 10.f);
  Filter ptfilter = aod::track::pt > cfgMinPt;
  Filter etafilter = aod::track::eta < cfgMaxEta;
  Filter properPIDfilter = aod::CMETrackPIDcolums::Npid != -1;

  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>>> Tracks_set1 = aod::CMETrackPIDcolums::Npid == 1;
  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>>> Tracks_set2 = aod::CMETrackPIDcolums::Npid == 2;
  Partition<soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>>> Tracks_set3 = aod::CMETrackPIDcolums::Npid == 3;
  void init(InitContext const&)
  {

    DetId = GetDetId(cfgDetName);
    RefAId = GetDetId(cfgRefAName);
    RefBId = GetDetId(cfgRefBName);

    if (DetId == RefAId || DetId == RefBId || RefAId == RefBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      DetId = 0;
      RefAId = 4;
      RefBId = 5;
    }

    AxisSpec axisCent{cfgaxisCent, "centrality"};
    AxisSpec axisQvec{cfgaxisQvec, "Q"};
    AxisSpec axisQvecF{cfgaxisQvecF, "Q"};
    AxisSpec axisEvtPl = {100, -1.0 * constants::math::PI, constants::math::PI};

    AxisSpec axisCos{cfgaxiscos, "angle function"};
    AxisSpec axisPt{cfgaxispt, "trasverse momentum"};
    AxisSpec axisCentMerged{cfgaxisCentMerged, "merged centrality"};

    AxisSpec axissumpt{cfgaxissumpt, "#it{p}_{T}^{sum}}"};
    AxisSpec axisdeltaeta{cfgaxisdeltaeta, "#Delta#eta"};
    AxisSpec axisvertexz = {100, -15., 15., "vrtx_{Z} [cm]"};

    histosQA.add(Form("QA/histEventCount"), "", {HistType::kTH1F, {{2, 0.0, 2.0}}});
    histosQA.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(1, "Not selected events");
    histosQA.get<TH1>(HIST("QA/histEventCount"))->GetXaxis()->SetBinLabel(2, "Selected events");
    histosQA.add(Form("QA/histVertexZRec"), "", {HistType::kTH1F, {axisvertexz}});
    histosQA.add(Form("QA/histCentrality"), "", {HistType::kTH1F, {axisCent}});
    histosQA.add(Form("QA/histQvec_CorrL0_V2"), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add(Form("QA/histQvec_CorrL1_V2"), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add(Form("QA/histQvec_CorrL2_V2"), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add(Form("QA/histQvec_CorrL3_V2"), "", {HistType::kTH3F, {axisQvecF, axisQvecF, axisCent}});
    histosQA.add(Form("QA/histEvtPl_CorrL0_V2"), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add(Form("QA/histEvtPl_CorrL1_V2"), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add(Form("QA/histEvtPl_CorrL2_V2"), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add(Form("QA/histEvtPl_CorrL3_V2"), "", {HistType::kTH2F, {axisEvtPl, axisCent}});
    histosQA.add(Form("QA/histQvecRes_SigRefAV2"), "", {HistType::kTH2F, {axisQvecF, axisCent}});
    histosQA.add(Form("QA/histQvecRes_SigRefBV2"), "", {HistType::kTH2F, {axisQvecF, axisCent}});
    histosQA.add(Form("QA/histQvecRes_RefARefBV2"), "", {HistType::kTH2F, {axisQvecF, axisCent}});

    histosQA.add(Form("V2/histCosDetV2"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("V2/histSinDetV2"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});

    histosQA.add(Form("V2/PID/histCosDetV2_Pi"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("V2/PID/histCosDetV2_Ka"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("V2/PID/histCosDetV2_Pr"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("V2/PID/histCosDetV2_Pi_Neg"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("V2/PID/histCosDetV2_Ka_Neg"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});
    histosQA.add(Form("V2/PID/histCosDetV2_Pr_Neg"), "", {HistType::kTH3F, {axisCentMerged, axisPt, axisCos}});

    if (OpenCME) {
      histosQA.add<TProfile>(Form("PIDCME/histgamama_PiKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
      histosQA.add<TProfile>(Form("PIDCME/histgamama_PiKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
      histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
      histosQA.add<TProfile>(Form("PIDCME/histgamama_PiPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
      histosQA.add<TProfile>(Form("PIDCME/histgamama_KaPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
      histosQA.add<TProfile>(Form("PIDCME/histgamama_KaPr_os"), "", {HistType::kTProfile, {axisCentMerged}});

      histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_ss"), "", {HistType::kTProfile, {axisCentMerged}});
      histosQA.add<TProfile>(Form("PIDCME/histdelta_PiKa_os"), "", {HistType::kTProfile, {axisCentMerged}});
      histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
      histosQA.add<TProfile>(Form("PIDCME/histdelta_PiPr_os"), "", {HistType::kTProfile, {axisCentMerged}});
      histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_ss"), "", {HistType::kTProfile, {axisCentMerged}});
      histosQA.add<TProfile>(Form("PIDCME/histdelta_KaPr_os"), "", {HistType::kTProfile, {axisCentMerged}});

      histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiKa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiKa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_PiPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_KaPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      histosQA.add<TProfile3D>(Form("PIDCME/Differential/histgamama_KaPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});

      histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiKa_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiKa_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_PiPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_KaPr_ss_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
      histosQA.add<TProfile3D>(Form("PIDCME/Differential/histdelta_KaPr_os_Dif"), "", {HistType::kTProfile3D, {axisCentMerged, axissumpt, axisdeltaeta}});
    }
  }

  template <typename CollType>
  bool SelEvent(const CollType& collision)
  {
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    if (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
    }
    return 1;
  }

  template <typename TrackType>
  bool SelTrack(const TrackType track)
  {
    if (!track.passedITSNCls())
      return false;
    if (!track.passedITSChi2NDF())
      return false;
    if (!track.passedITSHits())
      return false;
    if (!track.passedTPCCrossedRowsOverNCls())
      return false;
    if (!track.passedTPCChi2NDF())
      return false;
    if (!track.passedDCAxy())
      return false;
    if (!track.passedDCAz())
      return false;
    return true;
  }

  template <typename CollType>
  void fillHistosQvec(const CollType& collision, int nmode)
  {
    int DetInd = DetId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int RefAInd = RefAId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    int RefBInd = RefBId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    if (nmode == 2) {
      if (collision.qvecAmp()[DetId] > 1e-8) {
        histosQA.fill(HIST("QA/histQvec_CorrL0_V2"), collision.qvecRe()[DetInd], collision.qvecIm()[DetInd], collision.centFT0C());
        histosQA.fill(HIST("QA/histQvec_CorrL1_V2"), collision.qvecRe()[DetInd + 1], collision.qvecIm()[DetInd + 1], collision.centFT0C());
        histosQA.fill(HIST("QA/histQvec_CorrL2_V2"), collision.qvecRe()[DetInd + 2], collision.qvecIm()[DetInd + 2], collision.centFT0C());
        histosQA.fill(HIST("QA/histQvec_CorrL3_V2"), collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], collision.centFT0C());
        histosQA.fill(HIST("QA/histEvtPl_CorrL0_V2"), helperEP.GetEventPlane(collision.qvecRe()[DetInd], collision.qvecIm()[DetInd], nmode), collision.centFT0C());
        histosQA.fill(HIST("QA/histEvtPl_CorrL1_V2"), helperEP.GetEventPlane(collision.qvecRe()[DetInd + 1], collision.qvecIm()[DetInd + 1], nmode), collision.centFT0C());
        histosQA.fill(HIST("QA/histEvtPl_CorrL2_V2"), helperEP.GetEventPlane(collision.qvecRe()[DetInd + 2], collision.qvecIm()[DetInd + 2], nmode), collision.centFT0C());
        histosQA.fill(HIST("QA/histEvtPl_CorrL3_V2"), helperEP.GetEventPlane(collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], nmode), collision.centFT0C());
      }
      if (collision.qvecAmp()[DetId] > 1e-8 && collision.qvecAmp()[RefAId] > 1e-8 && collision.qvecAmp()[RefBId] > 1e-8) {
        histosQA.fill(HIST("QA/histQvecRes_SigRefAV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[RefAInd + 3], collision.qvecIm()[RefAInd + 3], nmode), nmode), collision.centFT0C());
        histosQA.fill(HIST("QA/histQvecRes_SigRefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[RefBInd + 3], collision.qvecIm()[RefBInd + 3], nmode), nmode), collision.centFT0C());
        histosQA.fill(HIST("QA/histQvecRes_RefARefBV2"), helperEP.GetResolution(helperEP.GetEventPlane(collision.qvecRe()[RefAInd + 3], collision.qvecIm()[RefAInd + 3], nmode), helperEP.GetEventPlane(collision.qvecRe()[RefBInd + 3], collision.qvecIm()[RefBInd + 3], nmode), nmode), collision.centFT0C());
      }
    }
  }

  template <typename CollType, typename TrackType>
  void fillHistosFlow_gamma_delta(const CollType& collision, const TrackType& track1, const TrackType& track2, const TrackType& track3, int nmode)
  {
    if (collision.qvecAmp()[DetId] < 1e-8) {
      return;
    }
    int DetInd = DetId * 4 + cfgnTotalSystem * 4 * (nmode - 2);
    bool kisPi = false, kisKa = false, kisPr = false;
    // bool kisPi_2 = false;
    bool kisKa_2 = false, kisPr_2 = false;
    float Psi_n = helperEP.GetEventPlane(collision.qvecRe()[DetInd + 3], collision.qvecIm()[DetInd + 3], nmode);
    for (auto& trk : track1) {
      if (!SelTrack(trk))
        continue;
      if (nmode == 2) {
        if (trk.sign() > 0) {
          histosQA.fill(HIST("V2/PID/histCosDetV2_Pi"), collision.centFT0C(), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() - Psi_n)));
        } else if (trk.sign() < 0) {
          histosQA.fill(HIST("V2/PID/histCosDetV2_Pi_Neg"), collision.centFT0C(), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() - Psi_n)));
        }
      }
    }
    for (auto& trk : track2) {
      if (!SelTrack(trk))
        continue;
      if (nmode == 2) {
        if (trk.sign() > 0) {
          histosQA.fill(HIST("V2/PID/histCosDetV2_Ka"), collision.centFT0C(), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() - Psi_n)));
        } else if (trk.sign() < 0) {
          histosQA.fill(HIST("V2/PID/histCosDetV2_Ka_Neg"), collision.centFT0C(), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() - Psi_n)));
        }
      }
    }
    for (auto& trk : track3) {
      if (!SelTrack(trk))
        continue;
      if (nmode == 2) {
        if (trk.sign() > 0) {
          histosQA.fill(HIST("V2/PID/histCosDetV2_Pr"), collision.centFT0C(), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() - Psi_n)));
        } else if (trk.sign() < 0) {
          histosQA.fill(HIST("V2/PID/histCosDetV2_Pr_Neg"), collision.centFT0C(), trk.pt(),
                        std::cos(static_cast<float>(nmode) * (trk.phi() - Psi_n)));
        }
      }
    }
    if (OpenCME) {
      for (auto& trk1 : track1) {
        for (auto& trk2 : track2) {
          if (trk1.globalIndex() == trk2.globalIndex())
            continue;
          if (nmode == 2) {
            if (trk1.sign() == trk2.sign()) {
              histosQA.fill(HIST("PIDCME/histgamama_PiKa_ss"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * Psi_n)));
              histosQA.fill(HIST("PIDCME/histdelta_PiKa_ss"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
              histosQA.fill(HIST("PIDCME/Differential/histgamama_PiKa_ss_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                            std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * Psi_n)));
              histosQA.fill(HIST("PIDCME/Differential/histdelta_PiKa_ss_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                            std::cos((trk1.phi() - trk2.phi())));
            } else {
              histosQA.fill(HIST("PIDCME/histgamama_PiKa_os"), collision.centFT0C(), std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * Psi_n)));
              histosQA.fill(HIST("PIDCME/histdelta_PiKa_os"), collision.centFT0C(), std::cos((trk1.phi() - trk2.phi())));
              histosQA.fill(HIST("PIDCME/Differential/histgamama_PiKa_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                            std::cos((trk1.phi() + trk2.phi() - static_cast<float>(nmode) * Psi_n)));
              histosQA.fill(HIST("PIDCME/Differential/histdelta_PiKa_os_Dif"), collision.centFT0C(), trk1.pt() + trk2.pt(), std::abs(trk1.eta() - trk2.eta()),
                            std::cos((trk1.phi() - trk2.phi())));
            }
          }
        }
      }
      for (auto& trk1 : track1) {
        for (auto& trk3 : track3) {
          if (trk1.globalIndex() == trk3.globalIndex())
            continue;
          if (nmode == 2) {
            if (trk1.sign() == trk3.sign()) {
              histosQA.fill(HIST("PIDCME/histgamama_PiPr_ss"), collision.centFT0C(), std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * Psi_n)));
              histosQA.fill(HIST("PIDCME/histdelta_PiPr_ss"), collision.centFT0C(), std::cos((trk1.phi() - trk3.phi())));
              histosQA.fill(HIST("PIDCME/Differential/histgamama_PiPr_ss_Dif"), collision.centFT0C(), trk1.pt() + trk3.pt(), std::abs(trk1.eta() - trk3.eta()),
                            std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * Psi_n)));
              histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPr_ss_Dif"), collision.centFT0C(), trk1.pt() + trk3.pt(), std::abs(trk1.eta() - trk3.eta()),
                            std::cos((trk1.phi() - trk3.phi())));
            } else {
              histosQA.fill(HIST("PIDCME/histgamama_PiPr_os"), collision.centFT0C(), std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * Psi_n)));
              histosQA.fill(HIST("PIDCME/histdelta_PiPr_os"), collision.centFT0C(), std::cos((trk1.phi() - trk3.phi())));
              histosQA.fill(HIST("PIDCME/Differential/histgamama_PiPr_os_Dif"), collision.centFT0C(), trk1.pt() + trk3.pt(), std::abs(trk1.eta() - trk3.eta()),
                            std::cos((trk1.phi() + trk3.phi() - static_cast<float>(nmode) * Psi_n)));
              histosQA.fill(HIST("PIDCME/Differential/histdelta_PiPr_os_Dif"), collision.centFT0C(), trk1.pt() + trk3.pt(), std::abs(trk1.eta() - trk3.eta()),
                            std::cos((trk1.phi() - trk3.phi())));
            }
          }
        }
      }
      for (auto& trk2 : track2) {
        for (auto& trk3 : track3) {
          if (trk2.globalIndex() == trk3.globalIndex())
            continue;
          if (nmode == 2) {
            if (trk2.sign() == trk3.sign()) {
              histosQA.fill(HIST("PIDCME/histgamama_KaPr_ss"), collision.centFT0C(), std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * Psi_n)));
              histosQA.fill(HIST("PIDCME/histdelta_KaPr_ss"), collision.centFT0C(), std::cos((trk2.phi() - trk3.phi())));
              histosQA.fill(HIST("PIDCME/Differential/histgamama_KaPr_ss_Dif"), collision.centFT0C(), trk2.pt() + trk3.pt(), std::abs(trk2.eta() - trk3.eta()),
                            std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * Psi_n)));
              histosQA.fill(HIST("PIDCME/Differential/histdelta_KaPr_ss_Dif"), collision.centFT0C(), trk2.pt() + trk3.pt(), std::abs(trk2.eta() - trk3.eta()),
                            std::cos((trk2.phi() - trk3.phi())));
            } else {
              histosQA.fill(HIST("PIDCME/histgamama_KaPr_os"), collision.centFT0C(), std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * Psi_n)));
              histosQA.fill(HIST("PIDCME/histdelta_KaPr_os"), collision.centFT0C(), std::cos((trk2.phi() - trk3.phi())));
              histosQA.fill(HIST("PIDCME/Differential/histgamama_KaPr_os_Dif"), collision.centFT0C(), trk2.pt() + trk3.pt(), std::abs(trk2.eta() - trk3.eta()),
                            std::cos((trk2.phi() + trk3.phi() - static_cast<float>(nmode) * Psi_n)));
              histosQA.fill(HIST("PIDCME/Differential/histdelta_KaPr_os_Dif"), collision.centFT0C(), trk2.pt() + trk3.pt(), std::abs(trk2.eta() - trk3.eta()),
                            std::cos((trk2.phi() - trk3.phi())));
            }
          }
        }
      }
    }
  }

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Qvectors>>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::Flags>> const& tracks)
  {
    histosQA.fill(HIST("QA/histEventCount"), 0.5);
    if (!SelEvent(collision)) {
      return;
    }
    histosQA.fill(HIST("QA/histEventCount"), 1.5);
    histosQA.fill(HIST("QA/histCentrality"), collision.centFT0C());
    histosQA.fill(HIST("QA/histVertexZRec"), collision.posZ());
    auto tracks1 = Tracks_set1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracks2 = Tracks_set2->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracks3 = Tracks_set3->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    mult1 = tracks1.size();
    mult2 = tracks2.size();
    mult3 = tracks3.size();
    if (mult1 < 1 || mult2 < 1 || mult3 < 1) // Reject Collisions without sufficient particles
      return;
    for (auto i = 0; i < static_cast<int>(cfgnMods->size()); i++) {
      int DetInd_global = DetId * 4 + cfgnTotalSystem * 4 * (cfgnMods->at(i) - 2);
      float Psi_n_global = helperEP.GetEventPlane(collision.qvecRe()[DetInd_global + 3], collision.qvecIm()[DetInd_global + 3], cfgnMods->at(i));
      for (auto& trk : tracks) {
        if (!SelTrack(trk))
          continue;
        histosQA.fill(HIST("V2/histSinDetV2"), collision.centFT0C(), trk.pt(),
                      std::sin(static_cast<float>(cfgnMods->at(i)) * (trk.phi() - Psi_n_global)));
        histosQA.fill(HIST("V2/histCosDetV2"), collision.centFT0C(), trk.pt(),
                      std::cos(static_cast<float>(cfgnMods->at(i)) * (trk.phi() - Psi_n_global)));
      }
      fillHistosQvec(collision, cfgnMods->at(i));
      fillHistosFlow_gamma_delta(collision, tracks1, tracks2, tracks3, cfgnMods->at(i));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FillPIDcolums>(cfgc),
    adaptAnalysisTask<pidcme>(cfgc),
  };
}
