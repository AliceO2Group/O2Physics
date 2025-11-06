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
// L* baryon v2
// Prottay Das (prottay.das@cern.ch)

#include "PWGLF/DataModel/EPCalibrationTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
struct lstarpbpbv2 {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

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
  ConfigurableAxis configresAxis{"configresAxis", {VARIABLE_WIDTH, -2.0, -1.5, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 1.5, 2.0}, "Resolution"};
  ConfigurableAxis configthnAxispT{"configthnAxisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configIMAxis{"configIMAxis", {VARIABLE_WIDTH, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.5, 1.7, 1.8, 1.9, 2.0}, "IM"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = nabs(aod::cent::centFT0C) < cfgCutCentrality;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>>;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec spAxis = {spNbins, lbinsp, hbinsp, "Sp"};

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{20, 0.0, 100.0}});

    histos.add("ResFT0CTPCSP", "ResFT0CTPCSP", HistType::kTHnSparseF, {configcentAxis, configresAxis});
    histos.add("ResFT0CFT0ASP", "ResFT0CFT0ASP", HistType::kTHnSparseF, {configcentAxis, configresAxis});
    histos.add("ResFT0ATPCSP", "ResFT0ATPCSP", HistType::kTHnSparseF, {configcentAxis, configresAxis});
    histos.add("hpv2vscentpt", "hpv2vscentpt", HistType::kTHnSparseF, {configIMAxis, configcentAxis, configthnAxispT, spAxis}, true);
    histos.add("hpv2vscentptrot", "hpv2vscentptrot", HistType::kTHnSparseF, {configIMAxis, configcentAxis, configthnAxispT, spAxis}, true);
    histos.add("hpv2vscentptlike", "hpv2vscentptlike", HistType::kTHnSparseF, {configIMAxis, configcentAxis, configthnAxispT, spAxis}, true);
  }

  double massKa = o2::constants::physics::MassKPlus;
  double massPr = o2::constants::physics::MassProton;

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
  bool strategySelectionPID(const T& candidate, int PID, int strategy)
  {
    if (PID == 0) {
      if (strategy == 0) {
        if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && candidate.hasTOF() && TMath::Abs(candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > 0.5) {
          return true;
        }
      } else if (strategy == 1) {
        if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && TMath::Sqrt(candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > 0.5) {
          return true;
        }
      } else if (strategy == 2) {
        if (candidate.pt() < 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && candidate.hasTOF() && TMath::Abs(candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > 0.5) {
          return true;
        }
        if (candidate.pt() >= 0.5 && TMath::Abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && !candidate.hasTOF()) {
          return true;
        }
      }
    }
    return false;
  }

  // TOF Veto
  template <typename T>
  bool selectionPID1(const T& candidate)
  {
    if (candidate.hasTOF()) {
      if (candidate.p() < 0.5 && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.p() >= 0.5 && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
        return true;
      }
    }
    if (!candidate.hasTOF()) {
      if (std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        return true;
      }
    }
    return false;
  }

  // TPC TOF
  template <typename T>
  bool selectionPID2(const T& candidate)
  {
    if (candidate.pt() < 0.7 && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.pt() >= 0.7 && candidate.pt() < 1.4) {
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
        return true;
      }
      if (!candidate.hasTOF()) {
        if (candidate.pt() >= 0.7 && candidate.pt() < 0.8 && candidate.tpcNSigmaPr() > -1.8 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.8 && candidate.pt() < 0.9 && candidate.tpcNSigmaPr() > -1.7 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.9 && candidate.pt() < 1.0 && candidate.tpcNSigmaPr() > -1.6 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 1.0 && candidate.pt() < 1.4 && candidate.tpcNSigmaPr() > -1.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
      }
    }
    if (candidate.pt() >= 1.4) {
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
        return true;
      }
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        return true;
      }
    }
    return false;
  }

  // TPC TOF
  template <typename T>
  bool selectionPID3(const T& candidate)
  {
    if (candidate.pt() < 0.7 && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.pt() >= 0.7 && candidate.pt() < 1.8) {
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
        return true;
      }
      if (!candidate.hasTOF()) {
        if (candidate.pt() >= 0.7 && candidate.pt() < 0.8 && candidate.tpcNSigmaPr() > -1.8 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.8 && candidate.pt() < 0.9 && candidate.tpcNSigmaPr() > -1.7 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.9 && candidate.pt() < 1.0 && candidate.tpcNSigmaPr() > -1.6 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
      }
    }
    if (candidate.pt() >= 1.8) {
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
        return true;
      }
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        return true;
      }
    }
    return false;
  }

  // TPC TOF
  template <typename T>
  bool selectionPID4(const T& candidate)
  {
    if (candidate.pt() < 0.7 && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.pt() >= 0.7 && candidate.pt() < 1.4) {
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
        return true;
      }
      if (!candidate.hasTOF()) {
        if (candidate.pt() >= 0.7 && candidate.pt() < 0.8 && candidate.tpcNSigmaPr() > -1.8 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.8 && candidate.pt() < 0.9 && candidate.tpcNSigmaPr() > -1.7 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.9 && candidate.pt() < 1.0 && candidate.tpcNSigmaPr() > -1.6 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 1.0 && candidate.pt() < 1.4 && candidate.tpcNSigmaPr() > -1.0 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
      }
    }
    if (candidate.pt() >= 1.4) {
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        if (candidate.pt() < 2.5 && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
          return true;
        }
        if (candidate.pt() >= 2.5 && candidate.pt() < 4 && candidate.tofNSigmaPr() > -2.0 && candidate.tofNSigmaPr() < nsigmaCutTOF) {
          return true;
        }
        if (candidate.pt() >= 4 && candidate.pt() < 5 && candidate.tofNSigmaPr() > -1.5 && candidate.tofNSigmaPr() < nsigmaCutTOF) {
          return true;
        }
        if (candidate.pt() >= 5 && candidate.tofNSigmaPr() > -1.0 && candidate.tofNSigmaPr() < nsigmaCutTOF) {
          return true;
        }
      }
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        return true;
      }
    }
    return false;
  }

  // TPC TOF
  template <typename T>
  bool selectionPID5(const T& candidate)
  {
    if (candidate.pt() < 0.7 && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.pt() >= 0.7 && candidate.pt() < 1.8) {
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
        return true;
      }
      if (!candidate.hasTOF()) {
        if (candidate.pt() >= 0.7 && candidate.pt() < 0.8 && candidate.tpcNSigmaPr() > -1.8 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.8 && candidate.pt() < 0.9 && candidate.tpcNSigmaPr() > -1.7 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.9 && candidate.pt() < 1.0 && candidate.tpcNSigmaPr() > -1.6 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
      }
    }
    if (candidate.pt() >= 1.8) {
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        if (candidate.pt() < 2.5 && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
          return true;
        }
        if (candidate.pt() >= 2.5 && candidate.pt() < 4 && candidate.tofNSigmaPr() > -2.0 && candidate.tofNSigmaPr() < nsigmaCutTOF) {
          return true;
        }
        if (candidate.pt() >= 4 && candidate.pt() < 5 && candidate.tofNSigmaPr() > -1.5 && candidate.tofNSigmaPr() < nsigmaCutTOF) {
          return true;
        }
        if (candidate.pt() >= 5 && candidate.tofNSigmaPr() > -1.0 && candidate.tofNSigmaPr() < nsigmaCutTOF) {
          return true;
        }
      }
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        return true;
      }
    }
    return false;
  }
  // TPC TOF
  template <typename T>
  bool selectionPID6(const T& candidate)
  {
    if (candidate.pt() < 0.7 && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.pt() >= 0.7) {
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
        if (candidate.pt() < 2.5 && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
          return true;
        }
        if (candidate.pt() >= 2.5 && candidate.pt() < 4 && candidate.tofNSigmaPr() > -2.0 && candidate.tofNSigmaPr() < nsigmaCutTOF) {
          return true;
        }
        if (candidate.pt() >= 4 && candidate.pt() < 5 && candidate.tofNSigmaPr() > -1.5 && candidate.tofNSigmaPr() < nsigmaCutTOF) {
          return true;
        }
        if (candidate.pt() >= 5 && candidate.tofNSigmaPr() > -1.0 && candidate.tofNSigmaPr() < nsigmaCutTOF) {
          return true;
        }
      }
      if (!candidate.hasTOF()) {
        if (candidate.pt() >= 0.7 && candidate.pt() < 0.8 && candidate.tpcNSigmaPr() > -1.8 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.8 && candidate.pt() < 0.9 && candidate.tpcNSigmaPr() > -1.7 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.9 && candidate.pt() < 1.0 && candidate.tpcNSigmaPr() > -1.6 && candidate.tpcNSigmaPr() < nsigmaCutTPC) {
          return true;
        }
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

  ROOT::Math::PxPyPzMVector LstarMother, daughter1, daughter2, kaonrot, lstarrot;

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!collision.sel8() || !collision.triggereventep() || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }

    o2::aod::ITSResponse itsResponse;

    auto centrality = collision.centFT0C();
    auto QFT0C = collision.qFT0C();
    auto QFT0A = collision.qFT0A();
    auto QTPC = collision.qTPC();
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();

    histos.fill(HIST("hCentrality"), centrality);

    histos.fill(HIST("ResFT0CTPCSP"), centrality, QFT0C * QTPC * TMath::Cos(2.0 * (psiFT0C - psiTPC)));
    histos.fill(HIST("ResFT0CFT0ASP"), centrality, QFT0C * QFT0A * TMath::Cos(2.0 * (psiFT0C - psiFT0A)));
    histos.fill(HIST("ResFT0ATPCSP"), centrality, QTPC * QFT0A * TMath::Cos(2.0 * (psiTPC - psiFT0A)));

    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      bool track1kaon = false;
      auto track1ID = track1.globalIndex();
      if (!strategySelectionPID(track1, 0, strategyPID)) {
        continue;
      }
      track1kaon = true;

      for (auto track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }

        if (track2.p() < 0.5 && !(itsResponse.nSigmaITS<o2::track::PID::Proton>(track2) > -nsigmaCutITS && itsResponse.nSigmaITS<o2::track::PID::Proton>(track2) < nsigmaCutITS)) { // required for protons only
          continue;
        }

        bool track2proton = false;
        if (PIDstrategy == 0 && !selectionPID1(track2)) {
          continue;
        }
        if (PIDstrategy == 1 && !selectionPID2(track2)) {
          continue;
        }
        if (PIDstrategy == 2 && !selectionPID3(track2)) {
          continue;
        }
        if (PIDstrategy == 3 && !selectionPID4(track2)) {
          continue;
        }
        if (PIDstrategy == 4 && !selectionPID5(track2)) {
          continue;
        }
        if (PIDstrategy == 5 && !selectionPID6(track2)) {
          continue;
        }

        track2proton = true;
        auto track2ID = track2.globalIndex();

        if (track2ID == track1ID) {
          continue;
        }
        if (!track1kaon || !track2proton) {
          continue;
        }
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPr);
        LstarMother = daughter1 + daughter2;
        if (TMath::Abs(LstarMother.Rapidity()) >= 0.5) {
          continue;
        }

        auto phiminuspsi = GetPhiInRange(LstarMother.Phi() - psiFT0C);
        auto v2 = TMath::Cos(2.0 * phiminuspsi) * QFT0C;

        // unlike sign
        if (track1.sign() * track2.sign() < 0) {
          histos.fill(HIST("hpv2vscentpt"), LstarMother.M(), centrality, LstarMother.Pt(), v2);
          if (fillRotation) {
            for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
              auto anglestart = confMinRot;
              auto angleend = confMaxRot;
              auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
              auto rotangle = anglestart + nrotbkg * anglestep;
              auto rotkaonPx = track1.px() * std::cos(rotangle) - track1.py() * std::sin(rotangle);
              auto rotkaonPy = track1.px() * std::sin(rotangle) + track1.py() * std::cos(rotangle);
              kaonrot = ROOT::Math::PxPyPzMVector(rotkaonPx, rotkaonPy, track1.pz(), massKa);
              lstarrot = kaonrot + daughter2;
              if (TMath::Abs(lstarrot.Rapidity()) >= 0.5) {
                continue;
              }

              auto phiminuspsirot = GetPhiInRange(lstarrot.Phi() - psiFT0C);
              auto v2rot = TMath::Cos(2.0 * phiminuspsirot) * QFT0C;

              histos.fill(HIST("hpv2vscentptrot"), lstarrot.M(), centrality, lstarrot.Pt(), v2rot);
            }
          }
        }
        // like sign
        if (track1.sign() * track2.sign() > 0) {
          histos.fill(HIST("hpv2vscentptlike"), lstarrot.M(), centrality, lstarrot.Pt(), v2);
        }
      }
    }
  }
  PROCESS_SWITCH(lstarpbpbv2, processSE, "Process Same event", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lstarpbpbv2>(cfgc, TaskName{"lstarpbpbv2"})};
}
