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
// Fast Lambda k0s eff QA task for correlation analysis
// prottay.das@cern.ch

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/SPCalibrationTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
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
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <deque>
#include <iostream>
#include <iterator>
#include <set> // <<< CHANGED: for dedup sets
#include <string>
#include <type_traits>
#include <unordered_map> // <<< CHANGED: for seenMap
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct lambdak0seff {

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  int mRunNumber;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;
  o2::ccdb::CcdbApi ccdbApi;
  TH1D* hwgtAL;
  // fill output
  struct : ConfigurableGroup {
    Configurable<bool> additionalEvSel{"additionalEvSel", false, "additionalEvSel"};
    Configurable<bool> additionalEvSel2{"additionalEvSel2", false, "additionalEvSel2"};
    Configurable<bool> additionalEvSel3{"additionalEvSel3", false, "additionalEvSel3"};
  } evselGrp;
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 50.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 30.0f, "Accepted minimum Centrality"};
  // proton track cut
  Configurable<float> cfgCutPT{"cfgCutPT", 0.15, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.1f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.1f, "DCAz range for tracks"};
  Configurable<int> cfgITScluster{"cfgITScluster", 5, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isPVContributor{"isPVContributor", true, "is PV contributor"};
  // Configs for V0
  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<float> ConfV0Rap{"ConfV0Rap", 0.8f, "Rapidity range of V0"};
  Configurable<double> ConfV0DCADaughMax{"ConfV0DCADaughMax", 0.2f, "Maximum DCA between the V0 daughters"};
  Configurable<double> ConfV0CPAMin{"ConfV0CPAMin", 0.9998f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 1.5f, "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 1.2, "Maximum V0 DCA to PV"};
  Configurable<double> cMinV0DCAPr{"cMinV0DCAPr", 0.05, "Minimum V0 daughters DCA to PV for Pr"};
  Configurable<double> cMinV0DCAPi{"cMinV0DCAPi", 0.05, "Minimum V0 daughters DCA to PV for Pi"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 20, "Maximum V0 life time"};
  Configurable<bool> analyzeLambda{"analyzeLambda", true, "flag for lambda analysis"};
  Configurable<bool> analyzeK0s{"analyzeK0s", false, "flag for K0s analysis"};
  Configurable<float> qtArmenterosMinForK0{"qtArmenterosMinForK0", 0.2, "Armenterous cut for K0s"};
  // config for V0 daughters
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.4, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2, "minimum daughter pion pt"};
  Configurable<float> rcrfc{"rcrfc", 0.8f, "Ratio of CR to FC"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 50.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 3, "PID selections for Lambda daughters"};

  struct : ConfigurableGroup {
    Configurable<int> IMNbins{"IMNbins", 100, "Number of bins in invariant mass"};
    Configurable<float> lbinIM{"lbinIM", 1.0, "lower bin value in IM histograms"};
    Configurable<float> hbinIM{"hbinIM", 1.2, "higher bin value in IM histograms"};
  } binGrp;
  struct : ConfigurableGroup {
    ConfigurableAxis configcentAxis{"configcentAxis", {VARIABLE_WIDTH, 0.0, 10.0, 40.0, 80.0, 150, 300}, "Cent FT0C"};
    ConfigurableAxis configthnAxispT{"configthnAxisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis configetaAxis{"configetaAxis", {VARIABLE_WIDTH, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8}, "Eta"};
    ConfigurableAxis configvzAxis{"configvzAxis", {VARIABLE_WIDTH, -10, -5, -0.0, 5, 10}, "Vz"};
  } axisGrp;

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec thnAxisInvMass{binGrp.IMNbins, binGrp.lbinIM, binGrp.hbinIM, "#it{M} (GeV/#it{c}^{2})"};

    std::vector<AxisSpec> runaxes2 = {thnAxisInvMass, axisGrp.configthnAxispT, axisGrp.configetaAxis, axisGrp.configvzAxis, axisGrp.configcentAxis};

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{axisGrp.configcentAxis}});
    histos.add("hCentralitymc", "Centrality distribution MC", kTH1F, {{axisGrp.configcentAxis}});
    histos.add("hSparseGenLambda", "hSparseGenLambda", HistType::kTHnSparseF, runaxes2, true);
    histos.add("hSparseGenAntiLambda", "hSparseGenAntiLambda", HistType::kTHnSparseF, runaxes2, true);
    histos.add("hSparseRecLambda", "hSparseRecLambda", HistType::kTHnSparseF, runaxes2, true);
    histos.add("hSparseRecAntiLambda", "hSparseRecAntiLambda", HistType::kTHnSparseF, runaxes2, true);
    histos.add("hSparseGenK0s", "hSparseGenK0s", HistType::kTHnSparseF, runaxes2, true);
    histos.add("hSparseRecK0s", "hSparseRecK0s", HistType::kTHnSparseF, runaxes2, true);

    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    LOGF(info, "Getting alignment offsets from the CCDB...");
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster && candidate.itsNClsInnerBarrel() >= 1)) {
      return false;
    }
    return true;
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate)
  {
    if (TMath::Abs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = TMath::Abs(candidate.dcaV0daughters());
    const float cpav0 = candidate.v0cosPA();

    float CtauLambda = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda;
    float CtauK0s = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massK0s;

    if (pT < ConfV0PtMin) {
      return false;
    }
    if (dcaDaughv0 > ConfV0DCADaughMax) {
      return false;
    }
    if (cpav0 < ConfV0CPAMin) {
      return false;
    }
    if (tranRad < ConfV0TranRadV0Min) {
      return false;
    }
    if (tranRad > ConfV0TranRadV0Max) {
      return false;
    }
    if (analyzeLambda && TMath::Abs(CtauLambda) > cMaxV0LifeTime) {
      return false;
    }
    if (analyzeK0s && TMath::Abs(CtauK0s) > cMaxV0LifeTime) {
      return false;
    }
    if (analyzeLambda && TMath::Abs(candidate.yLambda()) > ConfV0Rap) {
      return false;
    }
    if (analyzeK0s && TMath::Abs(candidate.yK0Short()) > ConfV0Rap) {
      return false;
    }
    return true;
  }

  template <typename V0, typename T>
  bool isSelectedV0Daughter(V0 const& candidate, T const& track, int pid, int pid2)
  {
    const auto tpcNClsF = track.tpcNClsFound();
    if (track.tpcNClsCrossedRows() < cfgTPCcluster) {
      return false;
    }

    if (tpcNClsF < ConfDaughTPCnclsMin) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < rcrfc) {
      return false;
    }

    if (analyzeLambda && pid == 0 && TMath::Abs(track.tpcNSigmaPr()) > ConfDaughPIDCuts) {
      return false;
    }
    if (analyzeLambda && pid == 1 && TMath::Abs(track.tpcNSigmaPi()) > ConfDaughPIDCuts) {
      return false;
    }
    if (analyzeK0s && TMath::Abs(track.tpcNSigmaPi()) > ConfDaughPIDCuts) {
      return false;
    }
    if (pid == 0 && (candidate.positivept() < cfgDaughPrPt || candidate.negativept() < cfgDaughPiPt)) {
      return false;
    }
    if (pid == 1 && (candidate.positivept() < cfgDaughPiPt || candidate.negativept() < cfgDaughPrPt)) {
      return false;
    }
    if (std::abs(candidate.positiveeta()) > ConfDaughEta || std::abs(candidate.negativeeta()) > ConfDaughEta) {
      return false;
    }

    if (analyzeLambda && pid2 == 0 && (TMath::Abs(candidate.dcapostopv()) < cMinV0DCAPr || TMath::Abs(candidate.dcanegtopv()) < cMinV0DCAPi)) {
      return false;
    }
    if (analyzeLambda && pid2 == 1 && (TMath::Abs(candidate.dcapostopv()) < cMinV0DCAPi || TMath::Abs(candidate.dcanegtopv()) < cMinV0DCAPr)) {
      return false;
    }
    if (analyzeK0s && (TMath::Abs(candidate.dcapostopv()) < cMinV0DCAPi || TMath::Abs(candidate.dcanegtopv()) < cMinV0DCAPi)) {
      return false;
    }
    if (analyzeK0s && (candidate.qtarm() / (std::abs(candidate.alpha()))) < 0.2) {
      return false;
    }

    return true;
  }

  bool shouldReject(bool LambdaTag, bool aLambdaTag,
                    const ROOT::Math::PxPyPzMVector& Lambdadummy,
                    const ROOT::Math::PxPyPzMVector& AntiLambdadummy)
  {
    const double minMass = 1.105;
    const double maxMass = 1.125;
    return (LambdaTag && aLambdaTag &&
            (Lambdadummy.M() > minMass && Lambdadummy.M() < maxMass) &&
            (AntiLambdadummy.M() > minMass && AntiLambdadummy.M() < maxMass));
  }

  ROOT::Math::PxPyPzMVector Lambda, AntiLambda, Lambdadummy, AntiLambdadummy, Proton, Pion, AntiProton, AntiPion, K0sdummy, K0s;
  double massLambda = o2::constants::physics::MassLambda;
  double massK0s = o2::constants::physics::MassK0Short;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);

  using CollisionMCTrueTable = aod::McCollisions;
  // using EventCandidatesMC = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::McCollisionLabels>>;
  using EventCandidatesMC = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
  using AllTrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa>>;
  using ResoV0s = aod::V0Datas;

  using TrackMCTrueTable = aod::McParticles;
  Preslice<ResoV0s> perCollision = aod::v0data::collisionId;

  ROOT::Math::PxPyPzMVector lambdadummymc, antiLambdadummymc, kshortdummymc, protonmc, pionmc, antiProtonmc, antiPionmc;

  void processMC(CollisionMCTrueTable::iterator const&, EventCandidatesMC const& RecCollisions, TrackMCTrueTable const& GenParticles, ResoV0s const& V0s, AllTrackCandidates const&)
  {

    for (auto& collision : RecCollisions) {

      if (!collision.sel8()) {
        continue;
      }
      double centrality = -999.;
      centrality = collision.centFT0C();
      double vz = collision.posZ();

      if (evselGrp.additionalEvSel && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
        continue;
      }
      if (evselGrp.additionalEvSel2 && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        continue;
      }
      if (evselGrp.additionalEvSel3 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
        continue;
      }

      histos.fill(HIST("hCentrality"), centrality);

      auto v0sThisColl = V0s.sliceBy(perCollision, collision.globalIndex());

      for (const auto& v0 : v0sThisColl) {

        auto postrack = v0.template posTrack_as<AllTrackCandidates>();
        auto negtrack = v0.template negTrack_as<AllTrackCandidates>();

        if (analyzeLambda && analyzeK0s)
          continue;
        if (!analyzeLambda && !analyzeK0s)
          continue;

        int LambdaTag = 0;
        int aLambdaTag = 0;
        int K0sTag = 0;

        const auto signpos = postrack.sign();
        const auto signneg = negtrack.sign();

        if (signpos < 0 || signneg > 0) {
          continue;
        }
        if (analyzeLambda) {
          if (isSelectedV0Daughter(v0, postrack, 0, 0) && isSelectedV0Daughter(v0, negtrack, 1, 0)) {
            LambdaTag = 1;
          }
          if (isSelectedV0Daughter(v0, negtrack, 0, 1) && isSelectedV0Daughter(v0, postrack, 1, 1)) {
            aLambdaTag = 1;
          }
        }
        if (analyzeK0s) {
          if (isSelectedV0Daughter(v0, postrack, 0, 0) && isSelectedV0Daughter(v0, negtrack, 1, 0)) {
            K0sTag = 1;
          }
        }

        if (analyzeLambda && (!LambdaTag && !aLambdaTag))
          continue;
        if (analyzeK0s && (!K0sTag))
          continue;

        if (!SelectionV0(collision, v0)) {
          continue;
        }

        if (analyzeLambda) {
          if (LambdaTag) {
            Proton = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPr);
            AntiPion = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPi);
            Lambdadummy = Proton + AntiPion;
          }
          if (aLambdaTag) {
            AntiProton = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPr);
            Pion = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPi);
            AntiLambdadummy = AntiProton + Pion;
          }

          if (shouldReject(LambdaTag, aLambdaTag, Lambdadummy, AntiLambdadummy)) {
            continue;
          }
        }

        if (analyzeK0s) {
          if (K0sTag) {
            Pion = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPi);
            AntiPion = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPi);
            K0sdummy = Pion + AntiPion;
          }
        }

        if (TMath::Abs(v0.eta()) > 0.8)
          continue;

        if (LambdaTag) {
          Lambda = Proton + AntiPion;
          histos.fill(HIST("hSparseRecLambda"), v0.mLambda(), v0.pt(), v0.eta(), vz, centrality);
        }
        if (aLambdaTag) {
          AntiLambda = AntiProton + Pion;
          histos.fill(HIST("hSparseRecAntiLambda"), v0.mAntiLambda(), v0.pt(), v0.eta(), vz, centrality);
        }
        if (K0sTag) {
          histos.fill(HIST("hSparseRecK0s"), v0.mK0Short(), v0.pt(), v0.eta(), vz, centrality);
        }
      }

      for (const auto& mcParticle : GenParticles) {

        if (analyzeLambda && std::abs(mcParticle.pdgCode()) != PDG_t::kLambda0) {
          continue;
        }
        if (analyzeK0s && std::abs(mcParticle.pdgCode()) != PDG_t::kK0Short) {
          continue;
        }
        if (std::abs(mcParticle.y()) > ConfV0Rap) {
          continue;
        }
        if (!mcParticle.isPhysicalPrimary() || !mcParticle.producedByGenerator()) {
          continue;
        }

        auto pdg1 = mcParticle.pdgCode();
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        int daughsize = 2;
        if (kDaughters.size() != daughsize) {
          continue;
        }

        int lambdacounter = 0;
        int antilambdacounter = 0;
        int k0scounter = 0;

        for (const auto& kCurrentDaughter : kDaughters) {
          if (std::abs(kCurrentDaughter.pdgCode()) != PDG_t::kProton && std::abs(kCurrentDaughter.pdgCode()) != PDG_t::kPiPlus) {
            continue;
          }
          if (kCurrentDaughter.pt() < 0.2 || TMath::Abs(kCurrentDaughter.eta()) > 0.8)
            continue;

          if (kCurrentDaughter.pdgCode() == PDG_t::kProton) {
            lambdacounter += 1;
            protonmc = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassProton);
          }
          if (kCurrentDaughter.pdgCode() == PDG_t::kPiMinus) {
            lambdacounter += 1;
            k0scounter += 1;
            antiPionmc = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassPionCharged);
          }

          if (kCurrentDaughter.pdgCode() == PDG_t::kProtonBar) {
            antilambdacounter += 1;
            antiProtonmc = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassProton);
          }
          if (kCurrentDaughter.pdgCode() == PDG_t::kPiPlus) {
            antilambdacounter += 1;
            k0scounter += 1;
            pionmc = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), o2::constants::physics::MassPionCharged);
          }
        }

        if (analyzeLambda && pdg1 == PDG_t::kLambda0 && lambdacounter == 2 && antilambdacounter != 2 && k0scounter != 2) {
          lambdadummymc = protonmc + antiPionmc;
          histos.fill(HIST("hSparseGenLambda"), lambdadummymc.M(), lambdadummymc.Pt(), lambdadummymc.Eta(), vz, centrality);
        }

        if (analyzeLambda && pdg1 == PDG_t::kLambda0Bar && antilambdacounter == 2 && lambdacounter != 2 && k0scounter != 2) {
          antiLambdadummymc = antiProtonmc + pionmc;
          histos.fill(HIST("hSparseGenAntiLambda"), antiLambdadummymc.M(), antiLambdadummymc.Pt(), antiLambdadummymc.Eta(), vz, centrality);
        }
        if (analyzeK0s && pdg1 == PDG_t::kK0Short && k0scounter == 2 && lambdacounter != 2 && antilambdacounter != 2) {
          kshortdummymc = antiPionmc + pionmc;
          histos.fill(HIST("hSparseGenK0s"), kshortdummymc.M(), kshortdummymc.Pt(), kshortdummymc.Eta(), vz, centrality);
        }
      }
    }
  }
  PROCESS_SWITCH(lambdak0seff, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdak0seff>(cfgc, TaskName{"lambdak0seff"})};
}
