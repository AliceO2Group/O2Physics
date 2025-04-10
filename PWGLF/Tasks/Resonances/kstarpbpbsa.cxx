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
// K* spin alignment
// Prottay Das (prottay.das@cern.ch)

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <string>
#include <vector>
#include <cmath>
#include <array>
#include <cstdlib>

#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include "TF1.h"

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct kstarpbpbsa {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Accepted maximum Centrality"};
  // track
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<float> nsigmaCutTOF{"nsigmacutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> removefaketrak{"removefaketrack", true, "Remove fake track from momentum difference"};
  Configurable<float> ConfFakeKaonCut{"ConfFakeKaonCut", 0.1, "Cut based on track from momentum difference"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<bool> OnlyTOF{"OnlyTOF", true, "OnlyTOF"};
  Configurable<int> strategyPID{"strategyPID", 2, "PID strategy"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0.0, "cut TOF beta"};
  Configurable<float> confMinRot{"confMinRot", 5.0 * TMath::Pi() / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * TMath::Pi() / 6.0, "Maximum of rotation"};
  Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};
  Configurable<bool> fillRotation{"fillRotation", true, "fill rotation"};
  Configurable<bool> like{"like", true, "fill rotation"};
  Configurable<int> spNbins{"spNbins", 2000, "Number of bins in sp"};
  Configurable<float> lbinsp{"lbinsp", -1.0, "lower bin value in sp histograms"};
  Configurable<float> hbinsp{"hbinsp", 1.0, "higher bin value in sp histograms"};
  Configurable<float> nsigmaCutITS{"nsigmaCutITS", 3.0, "Value of the ITS Nsigma cut"};
  Configurable<int> PIDstrategy{"PIDstrategy", 0, "0: TOF Veto, 1: TOF Veto opti, 2: TOF, 3: TOF loose 1, 4: TOF loose 2, 5: old pt dep"};

  ConfigurableAxis configcentAxis{"configcentAxis", {VARIABLE_WIDTH, 0.0, 10.0, 40.0, 80.0}, "Cent V0M"};
  ConfigurableAxis configrapAxis{"configrapAxis", {VARIABLE_WIDTH, -0.8, -0.4, 0.4, 0.8}, "Rapidity"};
  ConfigurableAxis configresAxis{"configresAxis", {VARIABLE_WIDTH, -2.0, -1.5, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 1.5, 2.0}, "Resolution"};
  ConfigurableAxis configthnAxispT{"configthnAxisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configIMAxis{"configIMAxis", {VARIABLE_WIDTH, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.5, 1.7, 1.8, 1.9, 2.0}, "IM"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = nabs(aod::cent::centFT0C) < cfgCutCentrality;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>>;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec spAxis = {spNbins, lbinsp, hbinsp, "Sp"};

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{20, 0.0, 100.0}});

    histos.add("ResFT0CTPC", "ResFT0CTPC", HistType::kTHnSparseF, {configcentAxis, configresAxis});
    histos.add("ResFT0CFT0A", "ResFT0CFT0A", HistType::kTHnSparseF, {configcentAxis, configresAxis});
    histos.add("ResFT0ATPC", "ResFT0ATPC", HistType::kTHnSparseF, {configcentAxis, configresAxis});
    histos.add("hSparseSAvsrapsameunlike", "hSparseSAvsrapsameunlike", HistType::kTHnSparseF, {configIMAxis, configthnAxispT, spAxis, configrapAxis, configcentAxis}, true);
    histos.add("hSparseSAvsrapsamelike", "hSparseSAvsrapsamelike", HistType::kTHnSparseF, {configIMAxis, configthnAxispT, spAxis, configrapAxis, configcentAxis}, true);
    histos.add("hSparseSAvsraprot", "hSparseSAvsraprot", HistType::kTHnSparseF, {configIMAxis, configthnAxispT, spAxis, configrapAxis, configcentAxis}, true);
  }

  double massKa = o2::constants::physics::MassKPlus;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPiMinus;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (useGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    if (!useGlobalTrack && !(candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool strategySelectionPID(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (candidate.pt() < 0.5 && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.pt() >= 0.5 && TMath::Sqrt(candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < nsigmaCutTOF) {
        return true;
      }
    } else if (PID == 1) {
      if (candidate.pt() < 0.5 && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.pt() >= 0.5 && TMath::Sqrt(candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi() + candidate.tofNSigmaPi() * candidate.tofNSigmaPi()) < nsigmaCutTOF) {
        return true;
      }
    }
    return false;
  }

  double GetPhiInRange(double phi)
  {
    double result = phi;
    while (result < 0) {
      result = result + 2. * TMath::Pi() / 2;
    }
    while (result > 2. * TMath::Pi() / 2) {
      result = result - 2. * TMath::Pi() / 2;
    }
    return result;
  }

  ROOT::Math::PxPyPzMVector KstarMother, daughter1, daughter2, fourVecDauCM;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY, eventplaneVec, eventplaneVecNorm, beamvector;
  ROOT::Math::PxPyPzMVector kstarrot, kaonrot, daughter2rot, fourVecDauCMrot;
  ROOT::Math::XYZVector threeVecDauCMrot, threeVecDauCMXYrot;

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!collision.sel8() || !collision.triggereventep() || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }

    auto centrality = collision.centFT0C();
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();

    histos.fill(HIST("hCentrality"), centrality);

    histos.fill(HIST("ResFT0CTPC"), centrality, TMath::Cos(2.0 * (psiFT0C - psiTPC)));
    histos.fill(HIST("ResFT0CFT0A"), centrality, TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
    histos.fill(HIST("ResFT0ATPC"), centrality, TMath::Cos(2.0 * (psiTPC - psiFT0A)));

    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      bool track1kaon = false;
      auto track1ID = track1.globalIndex();
      if (!strategySelectionPID(track1, 0)) {
        continue;
      }
      track1kaon = true;

      for (const auto& track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }

        bool track2pion = false;
        auto track2ID = track2.globalIndex();
        if (!strategySelectionPID(track2, 1)) {
          continue;
        }

        track2pion = true;
        if (track2ID == track1ID) {
          continue;
        }
        if (!track1kaon || !track2pion) {
          continue;
        }
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
        KstarMother = daughter1 + daughter2;
        if (std::abs(KstarMother.Rapidity()) > 0.9) {
          continue;
        }

        ROOT::Math::Boost boost{KstarMother.BoostToCM()};
        fourVecDauCM = boost(daughter1);
        threeVecDauCM = fourVecDauCM.Vect();
        threeVecDauCMXY = ROOT::Math::XYZVector(threeVecDauCM.X(), threeVecDauCM.Y(), 0.);
        eventplaneVec = ROOT::Math::XYZVector(std::cos(2.0 * psiFT0C), std::sin(2.0 * psiFT0C), 0);
        auto cosPhistarminuspsi = GetPhiInRange(fourVecDauCM.Phi() - psiFT0C);
        auto SA = TMath::Cos(2.0 * cosPhistarminuspsi);

        int track1Sign = track1.sign();
        int track2Sign = track2.sign();

        if (track1Sign * track2Sign < 0)
          histos.fill(HIST("hSparseSAvsrapsameunlike"), KstarMother.M(), KstarMother.Pt(), SA, KstarMother.Rapidity(), centrality);
        else if (track1Sign * track2Sign > 0)
          histos.fill(HIST("hSparseSAvsrapsamelike"), KstarMother.M(), KstarMother.Pt(), SA, KstarMother.Rapidity(), centrality);

        if (fillRotation) {
          for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
            auto anglestart = confMinRot;
            auto angleend = confMaxRot;
            auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
            auto rotangle = anglestart + nrotbkg * anglestep;
            auto rotkaonPx = track1.px() * std::cos(rotangle) - track1.py() * std::sin(rotangle);
            auto rotkaonPy = track1.px() * std::sin(rotangle) + track1.py() * std::cos(rotangle);
            kaonrot = ROOT::Math::PxPyPzMVector(rotkaonPx, rotkaonPy, track1.pz(), massKa);
            kstarrot = kaonrot + daughter2;

            if (std::abs(kstarrot.Rapidity()) > 0.9) {
              continue;
            }

            ROOT::Math::Boost boost{kstarrot.BoostToCM()};
            fourVecDauCMrot = boost(kaonrot);
            threeVecDauCMrot = fourVecDauCMrot.Vect();
            threeVecDauCMXYrot = ROOT::Math::XYZVector(threeVecDauCMrot.X(), threeVecDauCMrot.Y(), 0.);
            auto cosPhistarminuspsirot = GetPhiInRange(fourVecDauCMrot.Phi() - psiFT0C);
            auto SArot = TMath::Cos(2.0 * cosPhistarminuspsirot);

            histos.fill(HIST("hSparseSAvsraprot"), kstarrot.M(), kstarrot.Pt(), SArot, kstarrot.Rapidity(), centrality);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(kstarpbpbsa, processSE, "Process Same event", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kstarpbpbsa>(cfgc, TaskName{"kstarpbpbsa"})};
}
