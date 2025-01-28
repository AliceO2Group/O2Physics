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

/// \file flowPbpbPikp.cxx
/// \brief PID flow using the generic framework
/// \author Preet Bhanjan Pati <preet.bhanjan.pati@cern.ch>

#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include <vector>
#include <utility>
#include <array>
#include <string>

#include "Math/Vector4D.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/StepTHn.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Multiplicity.h"
#include "CommonConstants/PhysicsConstants.h"

#include "PWGCF/GenericFramework/Core/GFWPowerArray.h"
#include "PWGCF/GenericFramework/Core/GFW.h"
#include "PWGCF/GenericFramework/Core/GFWCumulant.h"
#include "PWGCF/GenericFramework/Core/FlowContainer.h"

#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"

#include <TProfile.h>
#include <TRandom3.h>
#include <TF1.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct FlowPbpbPikp {
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> noLaterThan{"noLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgTpcNsigmaCut, float, 3.0f, "TPC N-sigma cut for pions, kaons, protons")
  O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.5f, "Minimum pt to use TOF N-sigma")

  O2_DEFINE_CONFIGURABLE(cfgUseCosPA, bool, false, "Use Pointing angle for resonances")
  O2_DEFINE_CONFIGURABLE(cfgPhiCosPA, float, 0.04f, "Minimum Pointing angle for Phi")
  O2_DEFINE_CONFIGURABLE(cfgTpcCluster, int, 70, "Number of TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 2.0f, "DCAxy range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "DCAz range for tracks")
  O2_DEFINE_CONFIGURABLE(useGlobalTrack, bool, true, "use Global track")
  O2_DEFINE_CONFIGURABLE(cfgITScluster, int, 0, "Number of ITS cluster")
  O2_DEFINE_CONFIGURABLE(confRapidity, float, 0.5, "Rapidity cut")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancy, int, 3000, "Occupancy cut")
  O2_DEFINE_CONFIGURABLE(additionalEvsel, bool, false, "Additional event selcection")

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 5.00, 6.00, 8.00, 10.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisParticles{"axisParticles", {3, 0, 3}, "axis for different hadrons"};
  ConfigurableAxis axisPhiMass{"axisPhiMass", {50000, 0, 5}, "axis for invariant mass distibution for Phi"};
  ConfigurableAxis axisTPCsignal{"axisTPCsignal", {10000, 0, 1000}, "axis for TPC signal"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter =(nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz) && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;
  //using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  // using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::pidBayes, aod::pidBayesPi, aod::pidBayesKa, aod::pidBayesPr, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;

  SliceCache cache;
  Partition<AodTracks> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<AodTracks> negTracks = aod::track::signed1Pt < 0.0f;

  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TAxis* fPtAxis;
  TRandom3* fRndm = new TRandom3(0);

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(noLaterThan.value);

    histos.add("hPhi", "", {HistType::kTH1D, {axisPhi}});
    histos.add("hEta", "", {HistType::kTH1D, {axisEta}});
    histos.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
    histos.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    histos.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});
    histos.add("hPt", "", {HistType::kTH1D, {axisPt}});
    histos.add("c22_gap08", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c22_gap08_pi", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c22_gap08_ka", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c22_gap08_pr", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c24_full", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("KplusTPC", "", {HistType::kTH2D, {{axisPt, axisTPCsignal}}});
    histos.add("KminusTPC", "", {HistType::kTH2D, {{axisPt, axisTPCsignal}}});
    histos.add("TofTpcNsigma", "", {HistType::kTHnSparseD, {{axisParticles, axisNsigmaTPC, axisNsigmaTOF, axisPt}}});
    histos.add("partCount", "", {HistType::kTHnSparseD, {{axisParticles, axisMultiplicity, axisPt}}});
    histos.add("hPhiMass_sparse", "", {HistType::kTHnSparseD, {{axisPhiMass, axisPt, axisMultiplicity}}});

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    if (additionalEvsel) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x)", 0, 100);
      fMultCutLow->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x)", 0, 100);
      fMultCutHigh->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultMultPVCut = new TF1("fMultMultPVCut", "[0]+[1]*x+[2]*x*x", 0, 5000);
      fMultMultPVCut->SetParameters(-0.1, 0.785, -4.7e-05);
    }

    TObjArray* oba = new TObjArray();
    oba->Add(new TNamed("Ch08Gap22", "Ch08Gap22"));
    for (int i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Ch08Gap22_pt_%i", i + 1), "Ch08Gap22_pTDiff"));
    oba->Add(new TNamed("Pi08Gap22", "Pi08Gap22"));
    for (int i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Pi08Gap22_pt_%i", i + 1), "Pi08Gap22_pTDiff"));
    oba->Add(new TNamed("Ka08Gap22", "Ka08Gap22"));
    for (int i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Ka08Gap22_pt_%i", i + 1), "Ka08Gap22_pTDiff"));
    oba->Add(new TNamed("Pr08Gap22", "Pr08Gap22"));
    for (int i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Pr08Gap22_pt_%i", i + 1), "Pr08Gap22_pTDiff"));
    oba->Add(new TNamed("ChFull24", "ChFull24"));
    for (int i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("ChFull24_pt_%i", i + 1), "ChFull24_pTDiff"));

    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fPtAxis);
    fFC->Initialize(oba, axisMultiplicity, cfgNbootstrap);
    delete oba;

    fGFW->AddRegion("refN08", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP08", 0.4, 0.8, 1, 1);
    fGFW->AddRegion("full", -0.8, 0.8, 1, 512);
    fGFW->AddRegion("poi", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 1024);
    fGFW->AddRegion("ol", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 2048);

    // charged parts
    fGFW->AddRegion("poiN", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 128);
    fGFW->AddRegion("olN", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 256);

    // pion
    fGFW->AddRegion("poiNpi", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("olNpi", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 16);

    // kaon
    fGFW->AddRegion("poiNk", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 4);
    fGFW->AddRegion("olNk", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 32);

    // proton
    fGFW->AddRegion("poiNpr", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 8);
    fGFW->AddRegion("olNpr", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 64);

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Ch08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Pi08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Ka08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Pr08Gap22", kFALSE));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN refN08 | olN {2} refP08 {-2}", "Ch08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNpi refN08 | olNpi {2} refP08 {-2}", "Pi08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNk refN08 | olNk {2} refP08 {-2}", "Ka08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNpr refN08 | olNpr {2} refP08 {-2}", "Pr08Gap22", kTRUE));

    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poi full | ol {2 2 -2 -2}", "ChFull24", kTRUE));
    fGFW->CreateRegions();
  }

  enum Particles {
    PIONS,
    KAONS,
    PROTONS
  };

  template <typename TTrack>
  bool isFakeKaon(TTrack track)
  {
    const auto pglobal = track.p();
    const auto ptpc = track.tpcInnerParam();
    if (std::abs(pglobal - ptpc) > 0.1) {
      return true;
    }
    return false;
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const float& centrality)
  {
    auto multNTracksPV = collision.multNTracksPV();
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;

    return 1;
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (useGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTpcCluster)) {
      return false;
    }
    if (!useGlobalTrack && !(candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
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
    angle = std::acos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (cfgUseCosPA && angle < cfgPhiCosPA) {
      return false;
    }
    return true;
  }

  template <typename TTrack>
  int getNsigmaPID(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaCombined = {std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi()), std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa()), std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr())};
    int pid = -1;
    float nsigma = cfgTpcNsigmaCut;

    // Choose which nSigma to use
    std::array<float, 3> nSigmaToUse = (track.pt() > cfgTofPtCut && track.hasTOF()) ? nSigmaCombined : nSigmaTPC;

    // Select particle with the lowest nsigma
    for (int i = 0; i < 3; ++i) {
      if (std::abs(nSigmaToUse[i]) < nsigma) {
        pid = i;
        nsigma = std::abs(nSigmaToUse[i]);
      }
    }
    return pid + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
  }

  /*template <typename TTrack>
  std::pair<int, int> getBayesID(TTrack track)
  {
    std::array<int, 3> bayesprobs = {static_cast<int>(track.bayesPi()), static_cast<int>(track.bayesKa()), static_cast<int>(track.bayesPr())};
    int bayesid = -1;
    int prob = 0;

    for (int i = 0; i < 3; ++i) {
      if (bayesprobs[i] > prob && bayesprobs[i] > 80) {
        bayesid = i;
        prob = bayesprobs[i];
      }
    }
    return std::make_pair(bayesid, prob);
  }

  template <typename TTrack>
  int getBayesPIDIndex(TTrack track)
  {
    int maxProb[3] = {80, 80, 80};
    int pidID = -1;
    std::pair<int, int> idprob = getBayesID(track);
    if (idprob.first == PIONS || idprob.first == KAONS || idprob.first == PROTONS) { // 0 = pion, 1 = kaon, 2 = proton
      pidID = idprob.first;
      float nsigmaTPC[3] = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
      if (idprob.second > maxProb[pidID]) {
        if (std::fabs(nsigmaTPC[pidID]) > 3)
          return 0;
        return pidID + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
      } else {
        return 0;
      }
    }
    return 0;
  }*/

  template <typename TTrack, typename vector, char... chars>
  void resurrectParticle(TTrack trackplus, TTrack trackminus, vector plusdaug, vector minusdaug, vector mom, double plusmass, const ConstStr<chars...>& hist, const double cent)
  {
    for (auto const& [partplus, partminus] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(trackplus, trackminus))) {
      if (getNsigmaPID(partplus) != 2)
        continue;
      if (getNsigmaPID(partminus) != 2)
        continue;
      //if (partplus.globalIndex() == partminus.globalIndex())
      //  continue;
      if(isFakeKaon(partplus))
        continue;
      if(isFakeKaon(partminus))
        continue;
      if (!selectionPair(partplus, partminus))
        continue;
      if (!selectionTrack(partplus)) 
        continue;
      if (!selectionTrack(partminus)) 
        continue;
      
      plusdaug = ROOT::Math::PxPyPzMVector(partplus.px(), partplus.py(), partplus.pz(), plusmass);
      minusdaug = ROOT::Math::PxPyPzMVector(partminus.px(), partminus.py(), partminus.pz(), plusmass);
      mom = plusdaug + minusdaug;

      if(std::abs(mom.Rapidity()) < confRapidity) {
        histos.fill(hist, mom.M(), mom.Pt(), cent);
      }
    }
    return;
  }

  template <char... chars>
  void fillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0)
      return;
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        histos.fill(tarName, cent, val, dnx);
      return;
    }
    for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        histos.fill(tarName, fPtAxis->GetBinCenter(i), val, dnx);
    }
    return;
  }

  void fillFC(const GFW::CorrConfig& corrconf, const double& cent, const double& rndm)
  {
    double dnx, val;
    dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
    if (dnx == 0) {
      return;
    }
    if (!corrconf.pTDif) {
      val = fGFW->Calculate(corrconf, 0, kFALSE).real() / dnx;
      if (std::fabs(val) < 1) {
        fFC->FillProfile(corrconf.Head.c_str(), cent, val, dnx, rndm);
      }
      return;
    }
    for (int i = 1; i <= fPtAxis->GetNbins(); i++) {
      dnx = fGFW->Calculate(corrconf, i - 1, kTRUE).real();
      if (dnx == 0)
        continue;
      val = fGFW->Calculate(corrconf, i - 1, kFALSE).real() / dnx;
      if (std::fabs(val) < 1)
        fFC->FillProfile(Form("%s_pt_%i", corrconf.Head.c_str(), i), cent, val, dnx, rndm);
    }
    return;
  }

  ROOT::Math::PxPyPzMVector Phimom, kplusdaug, kminusdaug;
  double massKplus = o2::constants::physics::MassKPlus;

  void process(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks)
  {
    int nTot = tracks.size();
    if (nTot < 1)
      return;
    if (!collision.sel8() || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
      return;

    float lRandom = fRndm->Rndm();
    float vtxz = collision.posZ();
    const auto cent = collision.centFT0C();

    int occupancy = collision.trackOccupancyInTimeRange();

    if (occupancy > cfgCutOccupancy)
      return;

    if (additionalEvsel && !eventSelected(collision, cent))
      return;

    histos.fill(HIST("hVtxZ"), vtxz);
    histos.fill(HIST("hMult"), nTot);
    histos.fill(HIST("hCent"), collision.centFT0C());
    fGFW->Clear();
    

    auto posSlicedTracks = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negSlicedTracks = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    float weff = 1, wacc = 1;
    int pidIndex;

    std::cout<<" *************** Event check 7 *************"<<std::endl;

    resurrectParticle(posSlicedTracks, negSlicedTracks, kplusdaug, kminusdaug, Phimom, massKplus, HIST("hPhiMass_sparse"), cent);

    /*for (auto const& track1 : tracks) {
      double pt = track1.pt();
      histos.fill(HIST("hPhi"), track1.phi());
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hPt"), pt);

      histos.fill(HIST("TofTpcNsigma"), PIONS, track1.tpcNSigmaPi(), track1.tofNSigmaPi(), pt);
      histos.fill(HIST("TofTpcNsigma"), KAONS, track1.tpcNSigmaKa(), track1.tofNSigmaKa(), pt);
      histos.fill(HIST("TofTpcNsigma"), PROTONS, track1.tpcNSigmaPr(), track1.tofNSigmaPr(), pt);

      bool withinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
      bool withinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);       // within RF pT range

      // pidIndex = getBayesPIDIndex(track1);
      pidIndex = getNsigmaPID(track1);
      if (withinPtRef) {
        fGFW->Fill(track1.eta(), fPtAxis->FindBin(pt) - 1, track1.phi(), wacc * weff, 1);
        fGFW->Fill(track1.eta(), 1, track1.phi(), wacc * weff, 512);
      }
      if (withinPtPOI) {
        fGFW->Fill(track1.eta(), fPtAxis->FindBin(pt) - 1, track1.phi(), wacc * weff, 128);
        fGFW->Fill(track1.eta(), fPtAxis->FindBin(pt) - 1, track1.phi(), wacc * weff, 1024);
      }
      if (withinPtPOI && withinPtRef) {
        fGFW->Fill(track1.eta(), fPtAxis->FindBin(pt) - 1, track1.phi(), wacc * weff, 256);
        fGFW->Fill(track1.eta(), fPtAxis->FindBin(pt) - 1, track1.phi(), wacc * weff, 2048);
      }

      if (pidIndex) {
        histos.fill(HIST("partCount"), pidIndex - 1, cent, pt);
        if (withinPtPOI)
          fGFW->Fill(track1.eta(), fPtAxis->FindBin(pt) - 1, track1.phi(), wacc * weff, 1 << (pidIndex));
        if (withinPtPOI && withinPtRef)
          fGFW->Fill(track1.eta(), fPtAxis->FindBin(pt) - 1, track1.phi(), wacc * weff, 1 << (pidIndex + 3));
      }
    } // track1 loop ends

    // Filling c22 with ROOT TProfile
    fillProfile(corrconfigs.at(0), HIST("c22_gap08"), cent);
    fillProfile(corrconfigs.at(1), HIST("c22_gap08_pi"), cent);
    fillProfile(corrconfigs.at(2), HIST("c22_gap08_ka"), cent);
    fillProfile(corrconfigs.at(3), HIST("c22_gap08_pr"), cent);
    fillProfile(corrconfigs.at(4), HIST("c24_full"), cent);

    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      fillFC(corrconfigs.at(l_ind), cent, lRandom);
    }*/

  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FlowPbpbPikp>(cfgc)};
}
