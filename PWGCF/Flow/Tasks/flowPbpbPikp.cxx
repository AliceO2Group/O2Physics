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
#include <map>

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
#include "PWGCF/GenericFramework/Core/GFWWeights.h"
#include "PWGCF/GenericFramework/Core/GFWWeightsList.h"

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
  O2_DEFINE_CONFIGURABLE(cfgTpcCluster, int, 70, "Number of TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgFillWeights, bool, true, "Fill NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgOutputNUAWeights, bool, true, "Fill and output NUA weights")
  O2_DEFINE_CONFIGURABLE(cfgOutputRunByRun, bool, true, "Fill and output NUA weights run by run")
  O2_DEFINE_CONFIGURABLE(cfgEfficiency, std::string, "", "CCDB path to efficiency object")
  O2_DEFINE_CONFIGURABLE(cfgAcceptance, std::string, "", "CCDB path to acceptance object")
  O2_DEFINE_CONFIGURABLE(cfgTpcNsigmaCut, float, 2.0f, "TPC N-sigma cut for pions, kaons, protons")
  O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.5f, "Minimum pt to use TOF N-sigma")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 2.0f, "DCAxy range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "DCAz range for tracks")

  O2_DEFINE_CONFIGURABLE(cfgCutOccupancy, int, 3000, "Occupancy cut")
  O2_DEFINE_CONFIGURABLE(cfgUseGlobalTrack, bool, true, "use Global track")
  O2_DEFINE_CONFIGURABLE(cfgITScluster, int, 0, "Number of ITS cluster")
  O2_DEFINE_CONFIGURABLE(cfgTrackDensityCorrUse, bool, true, "Use track density efficiency correction")

  Configurable<std::vector<double>> cfgTrackDensityP0{"cfgTrackDensityP0", std::vector<double>{0.7217476707, 0.7384792571, 0.7542625668, 0.7640680200, 0.7701951667, 0.7755299053, 0.7805901710, 0.7849446786, 0.7957356586, 0.8113039262, 0.8211968966, 0.8280558878, 0.8329342135}, "parameter 0 for track density efficiency correction"};
  Configurable<std::vector<double>> cfgTrackDensityP1{"cfgTrackDensityP1", std::vector<double>{-2.169488e-05, -2.191913e-05, -2.295484e-05, -2.556538e-05, -2.754463e-05, -2.816832e-05, -2.846502e-05, -2.843857e-05, -2.705974e-05, -2.477018e-05, -2.321730e-05, -2.203315e-05, -2.109474e-05}, "parameter 1 for track density efficiency correction"};

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 5.00, 6.00, 8.00, 10.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisParticles{"axisParticles", {3, 0, 3}, "axis for different hadrons"};
  ConfigurableAxis axisTPCsignal{"axisTPCsignal", {10000, 0, 1000}, "axis for TPC signal"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz) && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;
  using AodTracksWithoutBayes = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;

  OutputObj<FlowContainer> fFC{FlowContainer("FlowContainer")};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  GFW* fGFW = new GFW();
  std::vector<GFW::CorrConfig> corrconfigs;
  TAxis* fPtAxis;
  TRandom3* fRndm = new TRandom3(0);

  std::map<int, std::vector<std::shared_ptr<TH3>>> th3sList;
  enum OutputSpecies {
    hRef = 0,
    hCharge,
    hPion,
    hKaon,
    hProton,
    kCount_OutputSpecies
  };
  int lastRunNumer = -1;
  std::vector<int> runNumbers;
  std::vector<GFWWeights*> mAcceptance;
  bool correctionsLoaded = false;

  // local track density correction
  std::vector<TF1*> funcEff;
  TH1D* hFindPtBin;
  TF1* funcV2;
  TF1* funcV3;
  TF1* funcV4;

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(noLaterThan.value);

    histos.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
    histos.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    histos.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});
    histos.add("hPhi", "", {HistType::kTH1D, {axisPhi}});
    histos.add("hPhiWeighted", "", {HistType::kTH1D, {axisPhi}});
    histos.add("hEta", "", {HistType::kTH1D, {axisEta}});
    histos.add("hPt", "", {HistType::kTH1D, {axisPt}});
    histos.add("c22_gap08", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c22_gap08_pi", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c22_gap08_ka", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c22_gap08_pr", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c24_full", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c24_gap08", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c24_gap08_pi", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c24_gap08_ka", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("c24_gap08_pr", "", {HistType::kTProfile, {axisMultiplicity}});
    histos.add("TofTpcNsigma", "", {HistType::kTHnSparseD, {{axisParticles, axisNsigmaTPC, axisNsigmaTOF, axisPt}}});
    histos.add("partCount", "", {HistType::kTHnSparseD, {{axisParticles, axisMultiplicity, axisPt}}});
    if (cfgOutputNUAWeights && !cfgOutputRunByRun) {
      histos.add<TH3>("NUA/hPhiEtaVtxz_ref", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      histos.add<TH3>("NUA/hPhiEtaVtxz_ch", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      histos.add<TH3>("NUA/hPhiEtaVtxz_pi", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      histos.add<TH3>("NUA/hPhiEtaVtxz_ka", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      histos.add<TH3>("NUA/hPhiEtaVtxz_pr", ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
    }

    o2::framework::AxisSpec axis = axisPt;
    int nPtBins = axis.binEdges.size() - 1;
    double* ptBins = &(axis.binEdges)[0];
    fPtAxis = new TAxis(nPtBins, ptBins);

    TObjArray* oba = new TObjArray();
    oba->Add(new TNamed("ChFull22", "ChFull22"));
    for (int i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("ChFull22_pt_%i", i + 1), "ChFull22_pTDiff"));

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

    oba->Add(new TNamed("Ch08Gap24", "Ch08Gap24"));
    for (int i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Ch08Gap24_pt_%i", i + 1), "Ch08Gap24_pTDiff"));
    oba->Add(new TNamed("Pi08Gap24", "Pi08Gap24"));
    for (int i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Pi08Gap24_pt_%i", i + 1), "Pi08Gap24_pTDiff"));
    oba->Add(new TNamed("Ka08Gap24", "Ka08Gap24"));
    for (int i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Ka08Gap24_pt_%i", i + 1), "Ka08Gap24_pTDiff"));
    oba->Add(new TNamed("Pr08Gap24", "Pr08Gap24"));
    for (int i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Pr08Gap24_pt_%i", i + 1), "Pr08Gap24_pTDiff"));

    oba->Add(new TNamed("Pi08BGap22", "Pi08BGap22"));
    for (int i = 0; i < fPtAxis->GetNbins(); i++)
      oba->Add(new TNamed(Form("Pi08BGap22_pt_%i", i + 1), "Pi08BGap22_pTDiff"));

    fFC->SetName("FlowContainer");
    fFC->SetXAxis(fPtAxis);
    fFC->Initialize(oba, axisMultiplicity, cfgNbootstrap);
    delete oba;

    // reference particles
    fGFW->AddRegion("refN08", -0.8, -0.4, 1, 1);
    fGFW->AddRegion("refP08", 0.4, 0.8, 1, 1);
    fGFW->AddRegion("full", -0.8, 0.8, 1, 512);

    // pt dependent charged particles
    fGFW->AddRegion("poiN", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 128);
    fGFW->AddRegion("olN", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 256);
    fGFW->AddRegion("poi", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 1024);
    fGFW->AddRegion("ol", -0.8, 0.8, 1 + fPtAxis->GetNbins(), 2048);

    // pion
    fGFW->AddRegion("poiNpi", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("olNpi", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 16);

    fGFW->AddRegion("poiPpi", 0.4, 0.8, 1 + fPtAxis->GetNbins(), 2);
    fGFW->AddRegion("olPpi", 0.4, 0.8, 1 + fPtAxis->GetNbins(), 16);

    // kaon
    fGFW->AddRegion("poiNk", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 4);
    fGFW->AddRegion("olNk", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 32);

    // proton
    fGFW->AddRegion("poiNpr", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 8);
    fGFW->AddRegion("olNpr", -0.8, -0.4, 1 + fPtAxis->GetNbins(), 64);

    // reference particles
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 -2}", "ChFull22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Ch08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Pi08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Ka08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Pr08Gap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("full {2 2 -2 -2}", "ChFull24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2 2} refP08 {-2 -2}", "Ch08Gap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2 2} refP08 {-2 -2}", "Pi08Gap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2 2} refP08 {-2 -2}", "Ka08Gap24", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2 2} refP08 {-2 -2}", "Pr08Gap24", kFALSE));

    // pt differential pois
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poi full | ol {2 -2}", "ChFull22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN refN08 | olN {2} refP08 {-2}", "Ch08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNpi refN08 | olNpi {2} refP08 {-2}", "Pi08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNk refN08 | olNk {2} refP08 {-2}", "Ka08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNpr refN08 | olNpr {2} refP08 {-2}", "Pr08Gap22", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poi full | ol {2 2 -2 -2}", "ChFull24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiN refN08 | olN {2 2} refP08 {-2 -2}", "Ch08Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNpi refN08 | olNpi {2 2} refP08 {-2 -2}", "Pi08Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNk refN08 | olNk {2 2} refP08 {-2 -2}", "Ka08Gap24", kTRUE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiNpr refN08 | olNpr {2 2} refP08 {-2 -2}", "Pr08Gap24", kTRUE));

    // Backward correlations
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("refN08 {2} refP08 {-2}", "Pi08BGap22", kFALSE));
    corrconfigs.push_back(fGFW->GetCorrelatorConfig("poiPpi refP08 | olPpi {2} refN08 {-2}", "Pi08BGap22", kTRUE));

    fGFW->CreateRegions();

    if (cfgTrackDensityCorrUse) {
      std::vector<double> pTEffBins = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0};
      hFindPtBin = new TH1D("hFindPtBin", "hFindPtBin", pTEffBins.size() - 1, &pTEffBins[0]);
      funcEff.resize(pTEffBins.size() - 1);
      // LHC24g3 Eff
      std::vector<double> f1p0 = cfgTrackDensityP0;
      std::vector<double> f1p1 = cfgTrackDensityP1;
      for (uint ifunc = 0; ifunc < pTEffBins.size() - 1; ifunc++) {
        funcEff[ifunc] = new TF1(Form("funcEff%i", ifunc), "[0]+[1]*x", 0, 3000);
        funcEff[ifunc]->SetParameters(f1p0[ifunc], f1p1[ifunc]);
      }
      funcV2 = new TF1("funcV2", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV2->SetParameters(0.0186111, 0.00351907, -4.38264e-05, 1.35383e-07, -3.96266e-10);
      funcV3 = new TF1("funcV3", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV3->SetParameters(0.0174056, 0.000703329, -1.45044e-05, 1.91991e-07, -1.62137e-09);
      funcV4 = new TF1("funcV4", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, 100);
      funcV4->SetParameters(0.008845, 0.000259668, -3.24435e-06, 4.54837e-08, -6.01825e-10);
    }
  }

  enum Particles {
    PIONS,
    KAONS,
    PROTONS
  };

  template <typename TTrack>
  bool selectionTrack(const TTrack& track)
  {
    if (cfgUseGlobalTrack && !(track.isGlobalTrack() && track.isPVContributor() && track.itsNCls() > cfgITScluster && track.tpcNClsFound() > cfgTpcCluster && track.hasTPC())) {
      return false;
    }
    if (!cfgUseGlobalTrack && !(track.isPVContributor() && track.itsNCls() > cfgITScluster && track.hasTPC())) {
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
    if (track.pt() >= cfgTofPtCut && !track.hasTOF())
      return -1;

    const int numSpecies = 3;
    // Select particle with the lowest nsigma
    for (int i = 0; i < numSpecies; ++i) {
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

    const int nParts = 3;
    const int bayesCut = 80;
    for (int i = 0; i < nParts; ++i) {
      if (bayesprobs[i] > prob && bayesprobs[i] > bayesCut) {
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
    const int nParts = 3;
    std::pair<int, int> idprob = getBayesID(track);
    if (idprob.first == PIONS || idprob.first == KAONS || idprob.first == PROTONS) { // 0 = pion, 1 = kaon, 2 = proton
      pidID = idprob.first;
      float nsigmaTPC[3] = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
      if (idprob.second > maxProb[pidID]) {
        if (std::fabs(nsigmaTPC[pidID]) > nParts)
          return 0;
        return pidID + 1; // shift the pid by 1, 1 = pion, 2 = kaon, 3 = proton
      } else {
        return 0;
      }
    }
    return 0;
  }*/

  template <char... chars>
  void fillProfile(const GFW::CorrConfig& corrconf, const ConstStr<chars...>& tarName, const double& cent)
  {
    double dnx, val;
    if (!corrconf.pTDif) {
      dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
      if (dnx == 0)
        return;
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
    if (!corrconf.pTDif) {
      dnx = fGFW->Calculate(corrconf, 0, kTRUE).real();
      if (dnx == 0) {
        return;
      }
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

  void createRunByRunHistos(int runNumber)
  {
    if (cfgOutputNUAWeights) {
      std::vector<std::shared_ptr<TH3>> tH3s(kCount_OutputSpecies);
      tH3s[hRef] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_ref", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      tH3s[hCharge] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_ch", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      tH3s[hPion] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_pi", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      tH3s[hKaon] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_ka", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      tH3s[hProton] = histos.add<TH3>(Form("NUA/%d/hPhiEtaVtxz_pr", runNumber), ";#varphi;#eta;v_{z}", {HistType::kTH3D, {axisPhi, {64, -1.6, 1.6}, {40, -10, 10}}});
      th3sList.insert(std::make_pair(runNumber, tH3s));
    }
  }

  void loadCorrections(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (correctionsLoaded)
      return;
    if (!cfgAcceptance.value.empty()) {
      uint64_t timestamp = bc.timestamp();
      mAcceptance.clear();
      mAcceptance.resize(kCount_OutputSpecies);
      mAcceptance[hRef] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_ref", timestamp);
      if (mAcceptance[hRef])
        LOGF(info, "Loaded acceptance weights from %s_ref (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hRef]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_ref (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hRef]);

      mAcceptance[hCharge] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_ch", timestamp);
      if (mAcceptance[hCharge])
        LOGF(info, "Loaded acceptance weights from %s_ch (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hCharge]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_ch (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hCharge]);

      mAcceptance[hPion] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_pi", timestamp);
      if (mAcceptance[hPion])
        LOGF(info, "Loaded acceptance weights from %s_pi (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hPion]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_pi (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hPion]);

      mAcceptance[hKaon] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_ka", timestamp);
      if (mAcceptance[hKaon])
        LOGF(info, "Loaded acceptance weights from %s_ka (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hKaon]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_ka (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hKaon]);

      mAcceptance[hProton] = ccdb->getForTimeStamp<GFWWeights>(cfgAcceptance.value + "_pr", timestamp);
      if (mAcceptance[hProton])
        LOGF(info, "Loaded acceptance weights from %s_pr (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hProton]);
      else
        LOGF(fatal, "Could not load acceptance weights from %s_pr (%p)", cfgAcceptance.value.c_str(), (void*)mAcceptance[hProton]);
    }

    correctionsLoaded = true;
  }

  template <typename TTrack>
  double getAcceptance(TTrack track, const double& vtxz, int index)
  { // 0 ref, 1 ch, 2 pi, 3 ka, 4 pr
    if (index < 0 || index >= kCount_OutputSpecies) {
      return 1;
    }
    double wacc = 1;
    if (!mAcceptance.empty() && correctionsLoaded) {
      if (!mAcceptance[index]) {
        LOGF(fatal, "Acceptance weights not loaded for index %d", index);
        return 1;
      }
      wacc = mAcceptance[index]->getNUA(track.phi(), track.eta(), vtxz);
    }
    return wacc;
  }

  template <typename TTrack>
  void fillWeights(const TTrack track, const double vtxz, const int& pid_index, const int& run)
  {
    double pt = track.pt();
    bool withinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
    bool withinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);       // within RF pT range

    if (cfgOutputRunByRun) {
      if (withinPtRef && !pid_index)
        th3sList[run][hRef]->Fill(track.phi(), track.eta(), vtxz); // pt-subset of charged particles for ref flow
      if (withinPtPOI)
        th3sList[run][hCharge + pid_index]->Fill(track.phi(), track.eta(), vtxz); // charged and id'ed particle weights
    } else {
      if (withinPtRef && !pid_index)
        histos.fill(HIST("NUA/hPhiEtaVtxz_ref"), track.phi(), track.eta(), vtxz); // pt-subset of charged particles for ref flow
      if (withinPtPOI) {
        switch (pid_index) {
          case 0:
            histos.fill(HIST("NUA/hPhiEtaVtxz_ch"), track.phi(), track.eta(), vtxz); // charged particle weights
            break;
          case 1:
            histos.fill(HIST("NUA/hPhiEtaVtxz_pi"), track.phi(), track.eta(), vtxz); // pion weights
            break;
          case 2:
            histos.fill(HIST("NUA/hPhiEtaVtxz_ka"), track.phi(), track.eta(), vtxz); // kaon weights
            break;
          case 3:
            histos.fill(HIST("NUA/hPhiEtaVtxz_pr"), track.phi(), track.eta(), vtxz); // proton weights
            break;
        }
      }
    }
  }

  void process(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracksWithoutBayes const& tracks)
  {
    int nTot = tracks.size();
    if (nTot < 1)
      return;

    if (!collision.sel8() || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))
      return;

    int occupancy = collision.trackOccupancyInTimeRange();
    if (occupancy > cfgCutOccupancy)
      return;

    float lRandom = fRndm->Rndm();
    float vtxz = collision.posZ();
    const auto cent = collision.centFT0C();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    int runNumber = bc.runNumber();
    if (cfgOutputRunByRun && runNumber != lastRunNumer) {
      lastRunNumer = runNumber;
      if (std::find(runNumbers.begin(), runNumbers.end(), runNumber) == runNumbers.end()) {
        // if run number is not in the preconfigured list, create new output histograms for this run
        createRunByRunHistos(runNumber);
        runNumbers.push_back(runNumber);
      }
    }

    histos.fill(HIST("hVtxZ"), vtxz);
    histos.fill(HIST("hMult"), nTot);
    histos.fill(HIST("hCent"), cent);
    fGFW->Clear();

    float weff = 1;
    int pidIndex;
    loadCorrections(bc); // load corrections for the each event

    // Track loop for calculating the Qn angles
    double psi2Est = 0, psi3Est = 0, psi4Est = 0;
    float wEPeff = 1;
    double v2 = 0, v3 = 0, v4 = 0;
    // be cautious, this only works for Pb-Pb
    // esimate the Qn angles and vn for this event
    if (cfgTrackDensityCorrUse) {
      double q2x = 0, q2y = 0;
      double q3x = 0, q3y = 0;
      double q4x = 0, q4y = 0;
      for (const auto& track : tracks) {
        bool withinPtRef = (cfgCutPtMin < track.pt()) && (track.pt() < cfgCutPtMax); // within RF pT range
        if (withinPtRef) {
          q2x += std::cos(2 * track.phi());
          q2y += std::sin(2 * track.phi());
          q3x += std::cos(3 * track.phi());
          q3y += std::sin(3 * track.phi());
          q4x += std::cos(4 * track.phi());
          q4y += std::sin(4 * track.phi());
        }
      }
      psi2Est = std::atan2(q2y, q2x) / 2.;
      psi3Est = std::atan2(q3y, q3x) / 3.;
      psi4Est = std::atan2(q4y, q4x) / 4.;

      v2 = funcV2->Eval(cent);
      v3 = funcV3->Eval(cent);
      v4 = funcV4->Eval(cent);
    }

    // Actual track loop
    for (auto const& track : tracks) {
      if (!selectionTrack(track))
        continue;
      double pt = track.pt();
      histos.fill(HIST("hPhi"), track.phi());
      histos.fill(HIST("hEta"), track.eta());
      histos.fill(HIST("hPt"), pt);

      histos.fill(HIST("TofTpcNsigma"), PIONS, track.tpcNSigmaPi(), track.tofNSigmaPi(), pt);
      histos.fill(HIST("TofTpcNsigma"), KAONS, track.tpcNSigmaKa(), track.tofNSigmaKa(), pt);
      histos.fill(HIST("TofTpcNsigma"), PROTONS, track.tpcNSigmaPr(), track.tofNSigmaPr(), pt);

      bool withinPtPOI = (cfgCutPtPOIMin < pt) && (pt < cfgCutPtPOIMax); // within POI pT range
      bool withinPtRef = (cfgCutPtMin < pt) && (pt < cfgCutPtMax);       // within RF pT range

      // pidIndex = getBayesPIDIndex(track);
      pidIndex = getNsigmaPID(track);

      weff = 1; // Initializing weff for each track
      // NUA weights
      if (cfgOutputNUAWeights)
        fillWeights(track, vtxz, pidIndex, runNumber);

      if (!withinPtPOI && !withinPtRef)
        return;
      double waccRef = getAcceptance(track, vtxz, 0);
      double waccPOI = withinPtPOI ? getAcceptance(track, vtxz, pidIndex + 1) : getAcceptance(track, vtxz, 0);
      if (withinPtRef && withinPtPOI && pidIndex)
        waccRef = waccPOI; // if particle is both (then it's overlap), override ref with POI

      // Track density correction
      if (cfgTrackDensityCorrUse && withinPtRef) {
        double fphi = v2 * std::cos(2 * (track.phi() - psi2Est)) + v3 * std::cos(3 * (track.phi() - psi3Est)) + v4 * std::cos(4 * (track.phi() - psi4Est));
        fphi = (1 + 2 * fphi);
        int pTBinForEff = hFindPtBin->FindBin(track.pt());
        if (pTBinForEff >= 1 && pTBinForEff <= hFindPtBin->GetNbinsX()) {
          wEPeff = funcEff[pTBinForEff - 1]->Eval(fphi * tracks.size());
          if (wEPeff > 0.) {
            wEPeff = 1. / wEPeff;
            weff *= wEPeff;
          }
        }
      } // end of track density correction loop

      if (withinPtRef) {
        histos.fill(HIST("hPhiWeighted"), track.phi(), waccRef);
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccRef * weff, 1);
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccRef * weff, 512);
      }
      if (withinPtPOI) {
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccPOI * weff, 128);
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccPOI * weff, 1024);
      }
      if (withinPtPOI && withinPtRef) {
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccPOI * weff, 256);
        fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccPOI * weff, 2048);
      }

      if (pidIndex) {
        histos.fill(HIST("partCount"), pidIndex - 1, cent, pt);
        if (withinPtPOI)
          fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccPOI * weff, 1 << (pidIndex));
        if (withinPtPOI && withinPtRef)
          fGFW->Fill(track.eta(), fPtAxis->FindBin(pt) - 1, track.phi(), waccPOI * weff, 1 << (pidIndex + 3));
      }
    } // track loop ends

    // Filling cumulants with ROOT TProfile
    fillProfile(corrconfigs.at(1), HIST("c22_gap08"), cent);
    fillProfile(corrconfigs.at(2), HIST("c22_gap08_pi"), cent);
    fillProfile(corrconfigs.at(3), HIST("c22_gap08_ka"), cent);
    fillProfile(corrconfigs.at(4), HIST("c22_gap08_pr"), cent);
    fillProfile(corrconfigs.at(5), HIST("c24_full"), cent);
    fillProfile(corrconfigs.at(6), HIST("c24_gap08"), cent);
    fillProfile(corrconfigs.at(7), HIST("c24_gap08_pi"), cent);
    fillProfile(corrconfigs.at(8), HIST("c24_gap08_ka"), cent);
    fillProfile(corrconfigs.at(9), HIST("c24_gap08_pr"), cent);

    for (uint l_ind = 0; l_ind < corrconfigs.size(); l_ind++) {
      fillFC(corrconfigs.at(l_ind), cent, lRandom);
    }

  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FlowPbpbPikp>(cfgc)};
}
