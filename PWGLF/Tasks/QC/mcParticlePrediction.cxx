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

///
/// \file   mcParticlePrediction.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Francesca Ercolessi francesca.ercolessi@cern.ch
/// \brief Task to build the predictions from the models based on the generated particles
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "PWGLF/Utils/mcParticle.h"
#include "PWGLF/Utils/inelGt.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "CommonConstants/LHCConstants.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pwglf;

// Particles
static const std::vector<std::string> parameterNames{"Enable"};
static constexpr int nParameters = 1;
static const int defaultParticles[PIDExtended::NIDsTot][nParameters]{{0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}};
bool enabledParticlesArray[PIDExtended::NIDsTot];

// Estimators
struct Estimators {
  typedef int estID;
  static constexpr estID FT0A = 0;
  static constexpr estID FT0C = 1;
  static constexpr estID FT0AC = 2;
  static constexpr estID FV0A = 3;
  static constexpr estID FDDA = 4;
  static constexpr estID FDDC = 5;
  static constexpr estID FDDAC = 6;
  static constexpr estID ZNA = 7;
  static constexpr estID ZNC = 8;
  static constexpr estID ZEM1 = 9;
  static constexpr estID ZEM2 = 10;
  static constexpr estID ZPA = 11;
  static constexpr estID ZPC = 12;
  static constexpr estID ITSIB = 13;
  static constexpr estID ETA05 = 14;
  static constexpr estID ETA08 = 15;
  static constexpr estID V0A = 16;  // (Run2)
  static constexpr estID V0C = 17;  // (Run2)
  static constexpr estID V0AC = 18; // (Run2 V0M)
  static constexpr estID nEstimators = 19;

  static constexpr const char* estimatorNames[nEstimators] = {"FT0A",
                                                              "FT0C",
                                                              "FT0AC",
                                                              "FV0A",
                                                              "FDDA",
                                                              "FDDC",
                                                              "FDDAC",
                                                              "ZNA",
                                                              "ZNC",
                                                              "ZEM1",
                                                              "ZEM2",
                                                              "ZPA",
                                                              "ZPC",
                                                              "ITSIB",
                                                              "ETA05",
                                                              "ETA08",
                                                              "V0A",
                                                              "V0C",
                                                              "V0AC"};
  static std::vector<std::string> arrayNames()
  {
    static std::vector<std::string> names;
    if (!names.empty()) {
      return names;
    }
    for (int i = 0; i < nEstimators; i++) {
      names.push_back(estimatorNames[i]);
    }
    return names;
  }
};
bool enabledEstimatorsArray[Estimators::nEstimators];
static const int defaultEstimators[Estimators::nEstimators][nParameters]{{0},  // FT0A
                                                                         {0},  // FT0C
                                                                         {1},  // FT0AC
                                                                         {0},  // FV0A
                                                                         {0},  // FDDA
                                                                         {0},  // FDDC
                                                                         {0},  // FDDAC
                                                                         {0},  // ZNA
                                                                         {0},  // ZNC
                                                                         {0},  // ZEM1
                                                                         {0},  // ZEM2
                                                                         {0},  // ZPA
                                                                         {0},  // ZPC
                                                                         {0},  // ITSIB
                                                                         {0},  // ETA05
                                                                         {0},  // ETA08
                                                                         {0},  // V0A (Run2)
                                                                         {0},  // V0C (Run2)
                                                                         {0}}; // V0AC (Run2 V0M)

// Histograms
std::array<std::shared_ptr<TH1>, Estimators::nEstimators> hestimators;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsVsITS;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsVsETA05;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsVsETA08;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvGenVsReco;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvGenVsReco_BCMC;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvGenVsRecoITS;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvRecoVsITS;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvRecoVsRecoITS;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvRecoVsRecoITS_BCMC;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvRecoVsFT0A;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvRecoVsBCId;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvVsBCId;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hvertexPosZ;
std::array<std::array<std::shared_ptr<TH2>, PIDExtended::NIDsTot>, Estimators::nEstimators> hpt;
std::array<std::array<std::shared_ptr<TH1>, PIDExtended::NIDsTot>, Estimators::nEstimators> hyield;

struct mcParticlePrediction {

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosRecoEvs{"HistosRecoEvs", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosYield{"HistosYield", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosPt{"HistosPt", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis binsEta{"binsEta", {100, -20, 20}, "Binning of the Eta axis"};
  ConfigurableAxis binsVxy{"binsVxy", {100, -10, 10}, "Binning of the production vertex (x and y) axis"};
  ConfigurableAxis binsVz{"binsVz", {100, -10, 10}, "Binning of the production vertex (z) axis"};
  ConfigurableAxis binsPt{"binsPt", {100, 0, 10}, "Binning of the Pt axis"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {300, -0.5, 299.5}, "Binning of the Multiplicity axis"};
  ConfigurableAxis binsMultiplicityReco{"binsMultiplicityReco", {1000, -0.5, -0.5 + 10000}, "Binning of the Multiplicity axis"};
  Configurable<LabeledArray<int>> enabledSpecies{"enabledSpecies",
                                                 {defaultParticles[0], PIDExtended::NIDsTot, nParameters, PIDExtended::arrayNames(), parameterNames},
                                                 "Particles enabled"};
  Configurable<LabeledArray<int>> enabledEstimators{"enabledEstimators",
                                                    {defaultEstimators[0], Estimators::nEstimators, nParameters, Estimators::arrayNames(), parameterNames},
                                                    "Estimators enabled"};
  Configurable<bool> selectInelGt0{"selectInelGt0", true, "Select only inelastic events"};
  Configurable<bool> selectPrimaries{"selectPrimaries", true, "Select only primary particles"};
  Configurable<bool> requireCoincidenceEstimators{"requireCoincidenceEstimators", false, "Asks for a coincidence when two estimators are used"};
  Configurable<bool> discardkIsGoodZvtxFT0vsPV{"discardkIsGoodZvtxFT0vsPV", false, "Select only collisions with matching BC and MC BC"};
  Configurable<bool> discardMismatchedBCs{"discardMismatchedBCs", false, "Select only collisions with matching BC and MC BC"};
  Configurable<bool> discardMismatchedFoundBCs{"discardMismatchedFoundBCs", false, "Select only collisions with matching found BC and MC BC"};
  Configurable<float> posZCut{"posZCut", 10.f, "Cut in the Z position of the primary vertex"};
  Configurable<float> collisionTimeResCut{"collisionTimeResCut", -40.f, "Cut in the collisionTimeRes"};
  Configurable<bool> requirekIsGoodZvtxFT0vsPV{"requirekIsGoodZvtxFT0vsPV", false, "Require kIsGoodZvtxFT0vsPV: small difference between z-vertex from PV and from FT0"};
  Configurable<bool> requirekIsVertexITSTPC{"requirekIsVertexITSTPC", false, "Require kIsVertexITSTPC: at least one ITS-TPC track (reject vertices built from ITS-only tracks)"};
  Configurable<bool> requirekIsVertexTOFmatched{"requirekIsVertexTOFmatched", false, "Require kIsVertexTOFmatched: at least one of vertex contributors is matched to TOF"};
  Configurable<bool> requirekIsVertexTRDmatched{"requirekIsVertexTRDmatched", false, "Require kIsVertexTRDmatched: at least one of vertex contributors is matched to TRD"};
  Configurable<bool> enableVsITSHistograms{"enableVsITSHistograms", true, "Enables the correlation between ITS and other estimators"};
  Configurable<bool> enableVsEta05Histograms{"enableVsEta05Histograms", true, "Enables the correlation between ETA05 and other estimators"};
  Configurable<bool> enableVsEta08Histograms{"enableVsEta08Histograms", true, "Enables the correlation between ETA08 and other estimators"};

  Service<o2::framework::O2DatabasePDG> pdgDB;
  o2::pwglf::ParticleCounter<o2::framework::O2DatabasePDG> mCounter;

  void init(o2::framework::InitContext&)
  {
    mCounter.mPdgDatabase = pdgDB.service;
    mCounter.mSelectPrimaries = selectPrimaries.value;
    const AxisSpec axisEta{binsEta, "#eta"};
    const AxisSpec axisVx{binsVxy, "Vx"};
    const AxisSpec axisVy{binsVxy, "Vy"};
    const AxisSpec axisVz{binsVz, "Vz"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisMultiplicity{binsMultiplicity, "Multiplicity (undefined)"};
    const AxisSpec axisMultiplicityReco{binsMultiplicityReco, "Multiplicity Reco. (undefined)"};
    const AxisSpec axisMultiplicityRecoITS{100, 0, 100, "Multiplicity Reco. ITSIB"};
    const AxisSpec axisMultiplicityGenV0s{100, 0, 100, "K0s gen"};
    const AxisSpec axisMultiplicityRecoV0s{20, 0, 20, "K0s reco"};
    const AxisSpec axisBCID{o2::constants::lhc::LHCMaxBunches, -0.5, -0.5 + o2::constants::lhc::LHCMaxBunches, "BC ID in orbit"};
    const AxisSpec axisBCIDMC{o2::constants::lhc::LHCMaxBunches, -0.5, -0.5 + o2::constants::lhc::LHCMaxBunches, "MC BC ID in orbit"};
    const AxisSpec axisFT0{1000, -5, 5, "Coll time FT0 (ps)"};

    auto h = histos.add<TH1>("collisions/generated", "collisions", kTH1D, {{10, -0.5, 9.5}});
    h->GetXaxis()->SetBinLabel(1, "Read");
    h->GetXaxis()->SetBinLabel(2, "INELgt0");
    h->GetXaxis()->SetBinLabel(3, "|Z|<10");
    h = histos.add<TH1>("collisions/reconstructed", "collisions", kTH1D, {{20, -0.5, 19.5}});
    h->GetXaxis()->SetBinLabel(1, "Read");
    h->GetXaxis()->SetBinLabel(2, "has_mcCollision");
    h->GetXaxis()->SetBinLabel(3, "sel8");
    h->GetXaxis()->SetBinLabel(4, "kIsBBT0A");
    h->GetXaxis()->SetBinLabel(5, "kIsBBT0C");
    h->GetXaxis()->SetBinLabel(6, "collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))");
    h->GetXaxis()->SetBinLabel(7, "globalBC == MC globalBC");
    h->GetXaxis()->SetBinLabel(8, "found globalBC == MC globalBC");
    h->GetXaxis()->SetBinLabel(9, "isINELgt0mc");
    h->GetXaxis()->SetBinLabel(10, "VTXz");
    h->GetXaxis()->SetBinLabel(11, "collisionTimeRes");
    h->GetXaxis()->SetBinLabel(11, "collisionTimeRes");
    h->GetXaxis()->SetBinLabel(12, "kIsGoodZvtxFT0vsPV");
    h->GetXaxis()->SetBinLabel(13, "kIsVertexITSTPC");
    h->GetXaxis()->SetBinLabel(14, "kIsVertexTOFmatched");
    h->GetXaxis()->SetBinLabel(15, "kIsVertexTRDmatched");

    histos.add("collisions/Reco/BCvsMCBC", "BC vs MC BC", kTH2D, {axisBCID, axisBCIDMC});
    histos.add<TH2>("collisions/Reco/FoundBCvsMCBC", "Found BC vs MC BC", kTH2D, {axisBCID, axisBCIDMC})->GetXaxis()->SetTitle("Found BC ID in orbit");
    histos.add<TH2>("collisions/Reco/FoundBCvsBC", "Found BC vs MC BC", kTH2D, {axisBCID, axisBCID})->GetXaxis()->SetTitle("Found BC ID in orbit");
    histos.add("collisions/Reco/collisionTime", "Collision Time", kTH1D, {{1000, -20, 20, "collisionTime"}});
    histos.add("collisions/Reco/collisionTimeRes", "Collision Time Res", kTH1D, {{1600, 0, 1600, "collisionTimeRes (ns)"}});
    histos.add("collisions/Reco/bcMinusfoundBc", "bcMinusfoundBc", kTH1D, {{1600, -1000, 1000, "bc - foundBc (ns)"}});
    histos.add("collisions/Reco/bcMinusfoundBcRatio", "bcMinusfoundBcRatio", kTH1D, {{1600, -40, 40, "(bc - foundBc)/collisionTimeRes"}});
    histos.add("collisions/Reco/bcMinusMcBcRatio", "bcMinusMcBcRatio", kTH1D, {{1600, -40, 40, "(bc - mcBc)/collisionTimeRes"}});
    histos.add("collisions/Reco/foundbcMinusMcBcRatio", "foundbcMinusMcBcRatio", kTH1D, {{1600, -40, 40, "(foundBc-mcBc)/collisionTimeRes"}});
    histos.add<TH1>("collisions/Reco/FT0A", "FT0A", kTH1D, {axisFT0})->GetXaxis()->SetTitle("Coll time FT0A (ps)");
    histos.add<TH1>("collisions/Reco/FT0C", "FT0C", kTH1D, {axisFT0})->GetXaxis()->SetTitle("Coll time FT0C (ps)");
    histos.add<TH1>("collisions/Reco/FT0AC", "FT0AC", kTH1D, {axisFT0})->GetXaxis()->SetTitle("Coll time FT0AC (ps)");
    histos.add("particles/eta/charged", "eta", kTH1D, {axisEta});
    histos.add("particles/eta/neutral", "eta", kTH1D, {axisEta});
    histos.add("particles/vtx/x", "Vx", kTH1D, {axisVx});
    histos.add("particles/vtx/y", "Vy", kTH1D, {axisVy});
    histos.add("particles/vtx/z", "Vz", kTH1D, {axisVz});
    histos.add("particles/FromCollVsFromMCColl", "FromCollVsFromMCColl", kTH2D, {{binsMultiplicity, "PV contributor particles (good bc)"}, {binsMultiplicityReco, "Particles in MC collision"}});
    histos.add("particles/FromCollVsFromMCCollBad", "FromCollVsFromMCCollBad", kTH2D, {{binsMultiplicity, "PV contributor particles (bad bc)"}, {binsMultiplicityReco, "Particles in MC collision"}});
    histos.add("particles/FromCollVsFromCollBad", "FromCollVsFromCollBad", kTH2D, {{binsMultiplicity, "PV contributor particles (good bc)"}, {binsMultiplicity, "PV contributor particles (bad bc)"}});
    histos.add("particles/FromCollBadOverFromCollVsVsFromMCColl", "FromCollBadOverFromCollVsVsFromMCColl", kTH2D, {{100, 0, 2, "bad/good"}, {binsMultiplicityReco, "Particles in MC collision"}});
    histos.add("V0s/V0RecovsPV", "V0s Reco + Ass vs PV", kTH2D, {axisMultiplicityRecoITS, axisMultiplicityRecoV0s});
    histos.add("V0s/V0RecoAssvsPV", "V0s Reco + Ass vs PV", kTH2D, {axisMultiplicityRecoITS, axisMultiplicityRecoV0s});
    histos.add("V0s/V0AssvsPV", "V0s Ass vs PV", kTH2D, {axisMultiplicityRecoITS, axisMultiplicityGenV0s});
    histos.add("V0s/V0RecoAssvsPV_TOFOneLeg", "V0s Reco + Ass + TOF 1 Leg vs PV", kTH2D, {axisMultiplicityRecoITS, axisMultiplicityRecoV0s});
    histos.add("V0s/V0RecoAssvsPV_TOFTwoLegs", "V0s Reco + Ass + TOF 2 Legs vs PV", kTH2D, {axisMultiplicityRecoITS, axisMultiplicityRecoV0s});

    for (int i = 0; i < Estimators::nEstimators; i++) {
      if (enabledEstimators->get(Estimators::estimatorNames[i], "Enable") != 1) {
        enabledEstimatorsArray[i] = false;
        continue;
      }
      LOG(info) << "Enabling estimator " << i << " " << Estimators::estimatorNames[i];
      enabledEstimatorsArray[i] = true;
    }

    h = histos.add<TH1>("particles/yields", "particles", kTH1D, {{PIDExtended::NIDsTot, -0.5, -0.5 + PIDExtended::NIDsTot}});
    for (int i = 0; i < PIDExtended::NIDsTot; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, PIDExtended::getName(i));
    }

    for (int i = 0; i < Estimators::nEstimators; i++) {
      if (!enabledEstimatorsArray[i]) {
        continue;
      }
      const char* name = Estimators::estimatorNames[i];
      hestimators[i] = histos.add<TH1>(Form("multiplicity/%s", name), name, kTH1D, {axisMultiplicity});
      hestimators[i]->GetXaxis()->SetTitle(Form("Multiplicity %s", name));

      auto make2DH = [&](const std::string& h, const char* ytitle) {
        auto hist = histos.add<TH2>(Form("%s%s", h.c_str(), name),
                                    name,
                                    kTH2D,
                                    {axisMultiplicity, axisMultiplicity});
        hist->GetXaxis()->SetTitle(Form("Multiplicity %s", name));
        hist->GetXaxis()->SetTitle(Form("Multiplicity %s", ytitle));
        return hist;
      };
      if (enableVsITSHistograms) {
        hestimatorsVsITS[i] = make2DH("multiplicity/vsITS/", Estimators::estimatorNames[Estimators::ITSIB]);
      }
      if (enableVsEta05Histograms) {
        hestimatorsVsETA05[i] = make2DH("multiplicity/vsETA05/", Estimators::estimatorNames[Estimators::ETA05]);
      }
      if (enableVsEta08Histograms) {
        hestimatorsVsETA08[i] = make2DH("multiplicity/vsETA08/", Estimators::estimatorNames[Estimators::ETA08]);
      }

      hvertexPosZ[i] = histos.add<TH2>(Form("multiplicity/posZ/%s", name), name, kTH2D, {{200, -20, 20, "pos Z"}, axisMultiplicity});
      hvertexPosZ[i]->GetYaxis()->SetTitle(Form("Multiplicity %s", name));

      if (!doprocessReco) { // Reco events
        continue;
      }

      hestimatorsRecoEvGenVsReco[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/GenVsReco/%s", name), name, kTH2D, {axisMultiplicity, axisMultiplicityReco});
      hestimatorsRecoEvGenVsReco[i]->GetXaxis()->SetTitle(Form("Multiplicity %s", name));
      hestimatorsRecoEvGenVsReco[i]->GetYaxis()->SetTitle(Form("Multiplicity Reco. %s", name));

      hestimatorsRecoEvGenVsReco_BCMC[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/GenVsReco_BCMC/%s", name), name, kTH2D, {axisMultiplicity, axisMultiplicityReco});
      hestimatorsRecoEvGenVsReco_BCMC[i]->GetXaxis()->SetTitle(Form("Multiplicity %s", name));
      hestimatorsRecoEvGenVsReco_BCMC[i]->GetYaxis()->SetTitle(Form("Multiplicity Reco. %s (BCMC)", name));

      hestimatorsRecoEvGenVsRecoITS[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/GenVsRecoITS/%s", name), name, kTH2D, {axisMultiplicity, axisMultiplicityRecoITS});
      hestimatorsRecoEvGenVsRecoITS[i]->GetXaxis()->SetTitle(Form("Multiplicity %s", name));

      hestimatorsRecoEvRecoVsITS[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/RecoVsITS/%s", name), name, kTH2D, {axisMultiplicityReco, axisMultiplicity});
      hestimatorsRecoEvRecoVsITS[i]->GetXaxis()->SetTitle(Form("Multiplicity Reco. %s", name));
      hestimatorsRecoEvRecoVsITS[i]->GetYaxis()->SetTitle(Form("Multiplicity %s", Estimators::estimatorNames[Estimators::ITSIB]));

      hestimatorsRecoEvRecoVsRecoITS[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/RecoVsRecoITS/%s", name), name, kTH2D, {axisMultiplicityReco, axisMultiplicityRecoITS});
      hestimatorsRecoEvRecoVsRecoITS[i]->GetXaxis()->SetTitle(Form("Multiplicity Reco. %s", name));

      hestimatorsRecoEvRecoVsRecoITS_BCMC[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/RecoVsRecoITS_BCMC/%s", name), name, kTH2D, {axisMultiplicityReco, axisMultiplicityRecoITS});
      hestimatorsRecoEvRecoVsRecoITS_BCMC[i]->GetXaxis()->SetTitle(Form("Multiplicity Reco. %s (BCMC)", name));

      hestimatorsRecoEvRecoVsFT0A[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/RecovsFT0A/%s", name), name, kTH2D, {axisMultiplicityReco, axisMultiplicity});
      hestimatorsRecoEvRecoVsFT0A[i]->GetXaxis()->SetTitle(Form("Multiplicity Reco. %s", name));
      hestimatorsRecoEvRecoVsFT0A[i]->GetYaxis()->SetTitle(Form("Multiplicity %s", Estimators::estimatorNames[Estimators::FT0A]));

      hestimatorsRecoEvRecoVsBCId[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/RecoVsBCId/%s", name), name, kTH2D, {axisBCID, axisMultiplicityReco});
      hestimatorsRecoEvRecoVsBCId[i]->GetYaxis()->SetTitle(Form("Multiplicity Reco. %s", name));

      hestimatorsRecoEvVsBCId[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/VsBCId/%s", name), name, kTH2D, {axisBCID, axisMultiplicity});
      hestimatorsRecoEvVsBCId[i]->GetYaxis()->SetTitle(Form("Multiplicity %s", name));
    }

    for (int i = 0; i < PIDExtended::NIDsTot; i++) {
      if (enabledSpecies->get(PIDExtended::getName(i), "Enable") != 1) {
        enabledParticlesArray[i] = false;
        continue;
      }
      LOG(info) << "Enabling particle " << i << " " << PIDExtended::getName(i);
      enabledParticlesArray[i] = true;
      for (int j = 0; j < Estimators::nEstimators; j++) {
        if (!enabledEstimatorsArray[j]) {
          continue;
        }
        const char* name = Estimators::estimatorNames[j];
        hpt[j][i] = histosPt.add<TH2>(Form("prediction/pt/%s/%s", name, PIDExtended::getName(i)), PIDExtended::getName(i), kTH2D, {axisPt, axisMultiplicity});
        hpt[j][i]->GetYaxis()->SetTitle(Form("Multiplicity %s", name));

        hyield[j][i] = histosYield.add<TH1>(Form("prediction/yield/%s/%s", name, PIDExtended::getName(i)), PIDExtended::getName(i), kTH1D, {axisMultiplicity});
        hyield[j][i]->GetYaxis()->SetTitle(Form("Multiplicity %s", name));
      }
    }
    histos.print();
    histosRecoEvs.print();
    histosPt.print();
    histosYield.print();
  }

  std::array<float, Estimators::nEstimators> genMult(const auto& mcParticles)
  {
    std::array<float, Estimators::nEstimators> nMult;
    if (enabledEstimatorsArray[Estimators::FT0A] || enabledEstimatorsArray[Estimators::FT0AC]) {
      nMult[Estimators::FT0A] = mCounter.countFT0A(mcParticles);
    }
    if (enabledEstimatorsArray[Estimators::FT0C] || enabledEstimatorsArray[Estimators::FT0AC]) {
      nMult[Estimators::FT0C] = mCounter.countFT0C(mcParticles);
    }
    if (enabledEstimatorsArray[Estimators::FT0AC]) {
      nMult[Estimators::FT0AC] = nMult[Estimators::FT0A] + nMult[Estimators::FT0C];
      if (requireCoincidenceEstimators && (nMult[Estimators::FT0A] <= 0.f || nMult[Estimators::FT0C] <= 0.f)) {
        nMult[Estimators::FT0AC] = 0;
      }
    }
    if (enabledEstimatorsArray[Estimators::FV0A]) {
      nMult[Estimators::FV0A] = mCounter.countFV0A(mcParticles);
    }
    if (enabledEstimatorsArray[Estimators::FDDA]) {
      nMult[Estimators::FDDA] = mCounter.countFDDA(mcParticles);
    }
    if (enabledEstimatorsArray[Estimators::FDDC]) {
      nMult[Estimators::FDDC] = mCounter.countFDDC(mcParticles);
    }
    if (enabledEstimatorsArray[Estimators::FDDAC]) {
      nMult[Estimators::FDDAC] = nMult[Estimators::FDDA] + nMult[Estimators::FDDC];
      if (requireCoincidenceEstimators && (nMult[Estimators::FDDA] <= 0.f || nMult[Estimators::FDDC] <= 0.f)) {
        nMult[Estimators::FDDAC] = 0;
      }
    }
    if (enabledEstimatorsArray[Estimators::ZNA]) {
      nMult[Estimators::ZNA] = mCounter.countZNA(mcParticles);
    }
    if (enabledEstimatorsArray[Estimators::ZNC]) {
      nMult[Estimators::ZNC] = mCounter.countZNC(mcParticles);
    }
    if (enabledEstimatorsArray[Estimators::ITSIB] || enableVsITSHistograms) {
      nMult[Estimators::ITSIB] = mCounter.countITSIB(mcParticles);
    }
    if (enabledEstimatorsArray[Estimators::ETA05] || enableVsEta05Histograms) {
      nMult[Estimators::ETA05] = mCounter.countEta05(mcParticles);
    }
    if (enabledEstimatorsArray[Estimators::ETA08] || enableVsEta08Histograms) {
      nMult[Estimators::ETA08] = mCounter.countEta08(mcParticles);
    }
    if (enabledEstimatorsArray[Estimators::V0A] || enabledEstimatorsArray[Estimators::V0AC]) {
      nMult[Estimators::V0A] = mCounter.countV0A(mcParticles);
    }
    if (enabledEstimatorsArray[Estimators::V0C] || enabledEstimatorsArray[Estimators::V0AC]) {
      nMult[Estimators::V0C] = mCounter.countV0C(mcParticles);
    }
    if (enabledEstimatorsArray[Estimators::V0AC]) {
      nMult[Estimators::V0AC] = nMult[Estimators::V0A] + nMult[Estimators::V0C];
      if (requireCoincidenceEstimators && (nMult[Estimators::V0A] <= 0 || nMult[Estimators::V0C] <= 0)) {
        nMult[Estimators::V0AC] = 0;
      }
    }
    return nMult;
  }

  void process(aod::McCollision const& mcCollision,
               aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("collisions/generated"), 0);
    if (selectInelGt0.value && !o2::pwglf::isINELgt0mc(mcParticles, pdgDB)) {
      return;
    }

    histos.fill(HIST("collisions/generated"), 1);
    if (std::abs(mcCollision.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("collisions/generated"), 2);

    const std::array<float, Estimators::nEstimators>& nMult = genMult(mcParticles);

    for (int i = 0; i < Estimators::nEstimators; i++) {
      if (!enabledEstimatorsArray[i]) {
        continue;
      }

      hestimators[i]->Fill(nMult[i]);
      if (enableVsITSHistograms) {
        hestimatorsVsITS[i]->Fill(nMult[i], nMult[Estimators::ITSIB]);
      }
      if (enableVsEta05Histograms) {
        hestimatorsVsETA05[i]->Fill(nMult[i], nMult[Estimators::ETA05]);
      }
      if (enableVsEta08Histograms) {
        hestimatorsVsETA08[i]->Fill(nMult[i], nMult[Estimators::ETA08]);
      }
      hvertexPosZ[i]->Fill(mcCollision.posZ(), nMult[i]);
    }

    for (const auto& particle : mcParticles) {
      particle.pdgCode();
      const auto id = PIDExtended::pdgToId(particle);
      if (id < 0) {
        continue;
      }
      if (!enabledParticlesArray[id]) {
        continue;
      }

      if (!particle.isPhysicalPrimary()) {
        continue;
      }

      TParticlePDG* p = pdgDB->GetParticle(particle.pdgCode());
      if (p) {
        if (std::abs(p->Charge()) > 1e-3) {
          histos.fill(HIST("particles/eta/charged"), particle.eta());
        } else {
          histos.fill(HIST("particles/eta/neutral"), particle.eta());
        }
      }

      if (std::abs(particle.y()) > 0.5) {
        continue;
      }

      histos.fill(HIST("particles/vtx/x"), particle.vx());
      histos.fill(HIST("particles/vtx/y"), particle.vy());
      histos.fill(HIST("particles/vtx/z"), particle.vz() - mcCollision.posZ());

      histos.fill(HIST("particles/yields"), id);
      for (int i = 0; i < Estimators::nEstimators; i++) {
        if (!enabledEstimatorsArray[i]) {
          continue;
        }
        hpt[i][id]->Fill(particle.pt(), nMult[i]);
        hyield[i][id]->Fill(nMult[i]);
      }
    }
  }

  using TracksMC = soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels>;

  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;
  SliceCache cache;
  void processReco(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::Mults, aod::EvSels, aod::FT0sCorrected>::iterator const& collision,
                   aod::McCollisions const& /*mcCollisions*/,
                   soa::Join<aod::BCs, aod::Run3MatchedToBCSparse> const& /*bcs*/,
                   aod::McParticles const& mcParticles,
                   TracksMC const& tracks,
                   aod::FT0s const&)
  {
    histos.fill(HIST("collisions/reconstructed"), 0);
    if (!collision.has_mcCollision()) {
      return;
    }
    const auto& mcCollision = collision.mcCollision();
    histos.fill(HIST("collisions/reconstructed"), 1);
    if (!collision.sel8()) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 2);
    if (!collision.selection_bit(aod::evsel::kIsBBT0A)) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 3);
    if (!collision.selection_bit(aod::evsel::kIsBBT0C)) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 4);
    if (discardkIsGoodZvtxFT0vsPV.value && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 5);

    const auto& recoBC = collision.bc_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>();
    const auto& foundBC = collision.foundBC_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>();
    const auto& mcBC = mcCollision.bc_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>();

    // Check that the BC in data and MC is the same
    if (discardMismatchedBCs.value && recoBC.globalBC() != mcBC.globalBC()) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 6);
    if (discardMismatchedFoundBCs.value && foundBC.globalBC() != mcBC.globalBC()) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 7);

    const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

    if (selectInelGt0.value && !o2::pwglf::isINELgt0mc(particlesInCollision, pdgDB)) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 8);

    if (std::abs(collision.posZ()) > posZCut.value) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 9);
    if (collisionTimeResCut.value > 0.f && collision.collisionTimeRes() > collisionTimeResCut.value) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 10);
    if (requirekIsGoodZvtxFT0vsPV.value && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 11);

    if (requirekIsVertexITSTPC.value && !collision.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 12);

    if (requirekIsVertexTOFmatched.value && !collision.selection_bit(aod::evsel::kIsVertexTOFmatched)) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 13);

    if (requirekIsVertexTRDmatched.value && !collision.selection_bit(aod::evsel::kIsVertexTRDmatched)) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 14);

    if (collision.t0ACorrectedValid()) {
      histos.fill(HIST("collisions/Reco/FT0A"), collision.t0ACorrected());
    }
    if (collision.t0CCorrectedValid()) {
      histos.fill(HIST("collisions/Reco/FT0C"), collision.t0CCorrected());
    }
    if (collision.t0ACValid()) {
      histos.fill(HIST("collisions/Reco/FT0AC"), collision.t0AC());
    }

    const auto& recoBCid = recoBC.globalBC() % o2::constants::lhc::LHCMaxBunches;
    const auto& mcBCid = mcBC.globalBC() % o2::constants::lhc::LHCMaxBunches;
    const auto& foundBCid = foundBC.globalBC() % o2::constants::lhc::LHCMaxBunches;
    const int diffRecoFoundBC = foundBC.globalBC() - recoBC.globalBC();
    const int diffRecoMCBC = recoBC.globalBC() - mcBC.globalBC();
    const int diffFoundMCBC = foundBC.globalBC() - mcBC.globalBC();
    histos.fill(HIST("collisions/Reco/BCvsMCBC"), recoBCid, mcBCid);
    histos.fill(HIST("collisions/Reco/FoundBCvsMCBC"), foundBCid, mcBCid);
    histos.fill(HIST("collisions/Reco/FoundBCvsBC"), foundBCid, recoBCid);
    histos.fill(HIST("collisions/Reco/bcMinusfoundBc"), (diffRecoFoundBC)*o2::constants::lhc::LHCBunchSpacingNS);
    histos.fill(HIST("collisions/Reco/bcMinusfoundBcRatio"), (diffRecoFoundBC)*o2::constants::lhc::LHCBunchSpacingNS / collision.collisionTimeRes());
    histos.fill(HIST("collisions/Reco/foundbcMinusMcBcRatio"), (diffFoundMCBC)*o2::constants::lhc::LHCBunchSpacingNS / collision.collisionTimeRes());
    histos.fill(HIST("collisions/Reco/bcMinusMcBcRatio"), (diffRecoMCBC)*o2::constants::lhc::LHCBunchSpacingNS / collision.collisionTimeRes());

    int particlesFromColl = 0;
    int particlesFromCollWrongBC = 0;
    for (const auto& track : tracks) {
      if (!track.isPVContributor()) {
        continue;
      }
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto& mcParticle = track.mcParticle();
      if (mcParticle.mcCollision().bc_as<soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>>().globalBC() == mcBC.globalBC()) {
        particlesFromColl++;
      } else {
        particlesFromCollWrongBC++;
      }
    }
    histos.fill(HIST("collisions/Reco/collisionTime"), collision.collisionTime());
    histos.fill(HIST("collisions/Reco/collisionTimeRes"), collision.collisionTimeRes());
    histos.fill(HIST("particles/FromCollVsFromMCColl"), particlesFromColl, particlesInCollision.size());
    histos.fill(HIST("particles/FromCollVsFromMCCollBad"), particlesFromCollWrongBC, particlesInCollision.size());
    histos.fill(HIST("particles/FromCollVsFromCollBad"), particlesFromColl, particlesFromCollWrongBC);
    histos.fill(HIST("particles/FromCollBadOverFromCollVsVsFromMCColl"), 1.f * particlesFromCollWrongBC / particlesFromColl, particlesInCollision.size());

    const std::array<float, Estimators::nEstimators>& nMult = genMult(particlesInCollision);

    float nMultReco[Estimators::nEstimators];
    nMultReco[Estimators::FT0A] = collision.multFT0A();
    nMultReco[Estimators::FT0C] = collision.multFT0C();
    nMultReco[Estimators::FT0AC] = collision.multFT0M();
    nMultReco[Estimators::FV0A] = collision.multFV0A();
    nMultReco[Estimators::FDDA] = collision.multFDDA();
    nMultReco[Estimators::FDDC] = collision.multFDDC();
    nMultReco[Estimators::FDDAC] = collision.multFDDM();
    nMultReco[Estimators::ZNA] = collision.multZNA();
    nMultReco[Estimators::ZNC] = collision.multZNC();
    nMultReco[Estimators::ITSIB] = collision.multNTracksPV();

    float nMultRecoMCBC[Estimators::nEstimators] = {0};
    if (mcBC.has_ft0()) {
      const auto& ft0 = mcBC.ft0();
      for (auto amplitude : ft0.amplitudeA()) {
        nMultRecoMCBC[Estimators::FT0A] += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        nMultRecoMCBC[Estimators::FT0C] += amplitude;
      }
      nMultRecoMCBC[Estimators::FT0AC] = nMultRecoMCBC[Estimators::FT0A] + nMultRecoMCBC[Estimators::FT0C];
    } else {
      nMultRecoMCBC[Estimators::FT0A] = -999.f;
      nMultRecoMCBC[Estimators::FT0C] = -999.f;
    }

    for (int i = 0; i < Estimators::nEstimators; i++) {
      if (!enabledEstimatorsArray[i]) {
        continue;
      }
      hestimatorsRecoEvGenVsReco[i]->Fill(nMult[i], nMultReco[i]);
      hestimatorsRecoEvGenVsReco_BCMC[i]->Fill(nMult[i], nMultRecoMCBC[i]);
      hestimatorsRecoEvGenVsRecoITS[i]->Fill(nMult[i], nMultReco[Estimators::ITSIB]);
      hestimatorsRecoEvRecoVsITS[i]->Fill(nMultReco[i], nMult[Estimators::ITSIB]);
      hestimatorsRecoEvRecoVsRecoITS[i]->Fill(nMultReco[i], nMultReco[Estimators::ITSIB]);
      hestimatorsRecoEvRecoVsRecoITS_BCMC[i]->Fill(nMultRecoMCBC[i], nMultReco[Estimators::ITSIB]);
      hestimatorsRecoEvRecoVsFT0A[i]->Fill(nMultReco[i], nMult[Estimators::FT0A]);
      hestimatorsRecoEvRecoVsBCId[i]->Fill(foundBCid, nMult[i]);
      hestimatorsRecoEvVsBCId[i]->Fill(foundBCid, nMultReco[i]);
    }
  }
  PROCESS_SWITCH(mcParticlePrediction, processReco, "Process the reco info", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<mcParticlePrediction>(cfgc)}; }
