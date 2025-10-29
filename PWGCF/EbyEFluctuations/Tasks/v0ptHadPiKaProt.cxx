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

/// \file v0ptHadPiKaProt.cxx
/// \brief Task for analyzing v0(pT) of inclusive hadrons, pions, kaons, and, protons
/// \author Swati Saha

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"
#include <CCDB/BasicCCDBManager.h>

#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TH2F.h>
#include <THn.h>
#include <TList.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static constexpr float LongArrayFloat[3][20] = {{1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {2.1, 2.2, 2.3, -2.1, -2.2, -2.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}, {3.1, 3.2, 3.3, -3.1, -3.2, -3.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2, 1.3, -1.1, -1.2, -1.3, 1.1, 1.2}};

struct V0ptHadPiKaProt {

  // ITS response
  o2::aod::ITSResponse itsResponse;
  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  enum Particles {
    PIONS = 0,
    KAONS,
    PROTONS
  };
  enum ParticleNsigma {
    kPionUpCut = 0,
    kKaonUpCut,
    kProtonUpCut,
    kPionLowCut,
    kKaonLowCut,
    kProtonLowCut
  };
  enum DetectorType {
    kTPC = 0,
    kTOF,
    kITS
  };
  enum CentralityEstimator {
    kFT0C = 0,
    kFT0A,
    kFT0M,
    kFV0A
  };

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutTpcChi2NCl{"cfgCutTpcChi2NCl", 2.5f, "Maximum TPCchi2NCl"};
  Configurable<float> cfgCutItsChi2NCl{"cfgCutItsChi2NCl", 36.0f, "Maximum ITSchi2NCl"};
  Configurable<float> cfgCutTrackDcaZ{"cfgCutTrackDcaZ", 2.0f, "Maximum DcaZ"};
  Configurable<int> cfgITScluster{"cfgITScluster", 1, "Minimum Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 80, "Minimum Number of TPC cluster"};
  Configurable<int> cfgTPCnCrossedRows{"cfgTPCnCrossedRows", 70, "Minimum Number of TPC crossed-rows"};
  Configurable<float> cfgCutPtUpperTPC{"cfgCutPtUpperTPC", 0.6f, "Upper pT cut for PID using TPC only"};
  Configurable<float> cfgnSigmaOtherParticles{"cfgnSigmaOtherParticles", 3.0f, "PID nSigma cut to remove other particles (default:3)"};
  Configurable<float> cfgnSigmaCutTPC{"cfgnSigmaCutTPC", 2.0f, "PID nSigma cut for TPC"};
  Configurable<float> cfgnSigmaCutTOF{"cfgnSigmaCutTOF", 2.0f, "PID nSigma cut for TOF"};
  Configurable<float> cfgnSigmaCutCombTPCTOF{"cfgnSigmaCutCombTPCTOF", 2.0f, "PID nSigma combined cut for TPC and TOF"};
  ConfigurableAxis nchAxis{"nchAxis", {5000, 0.5, 5000.5}, ""};
  ConfigurableAxis centAxis{"centAxis", {90, 0., 90.}, "Centrality/Multiplicity percentile bining"};
  Configurable<float> cfgCutPtLower{"cfgCutPtLower", 0.2f, "Lower pT cut"};
  Configurable<float> cfgCutPtLowerProt{"cfgCutPtLowerProt", 0.2f, "Lower pT cut"};
  Configurable<float> cfgCutPtUpper{"cfgCutPtUpper", 10.0f, "Higher pT cut for inclusive hadron analysis"};
  Configurable<float> cfgCutPtUpperPID{"cfgCutPtUpperPID", 6.0f, "Higher pT cut for identified particle analysis"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "absolute Eta cut"};
  Configurable<float> cfgCutEtaLeft{"cfgCutEtaLeft", 0.8f, "Left end of eta gap"};
  Configurable<float> cfgCutEtaRight{"cfgCutEtaRight", 0.8f, "Right end of eta gap"};
  Configurable<int> cfgNSubsample{"cfgNSubsample", 10, "Number of subsamples"};
  Configurable<int> cfgCentralityChoice{"cfgCentralityChoice", 0, "Which centrality estimator? 0-->FT0C, 1-->FT0A, 2-->FT0M, 3-->FV0A"};
  Configurable<bool> cfgEvSelkNoSameBunchPileup{"cfgEvSelkNoSameBunchPileup", true, "Pileup removal"};
  Configurable<bool> cfgUseGoodITSLayerAllCut{"cfgUseGoodITSLayerAllCut", true, "Remove time interval with dead ITS zone"};
  Configurable<bool> cfgEvSelkNoITSROFrameBorder{"cfgEvSelkNoITSROFrameBorder", true, "ITSROFrame border event selection cut"};
  Configurable<bool> cfgEvSelkNoTimeFrameBorder{"cfgEvSelkNoTimeFrameBorder", true, "TimeFrame border event selection cut"};
  Configurable<bool> cfgEvSelUseGoodZvtxFT0vsPV{"cfgEvSelUseGoodZvtxFT0vsPV", true, "GoodZvertex and FT0 vs PV cut"};
  Configurable<bool> cfgUseItsPID{"cfgUseItsPID", false, "Use ITS PID for particle identification"};
  Configurable<float> cfgPtCutTOF{"cfgPtCutTOF", 0.3f, "Minimum pt to use TOF N-sigma"};
  Configurable<LabeledArray<float>> nSigmas{"nSigmas", {LongArrayFloat[0], 3, 6, {"TPC", "TOF", "ITS"}, {"pos_pi", "pos_ka", "pos_pr", "neg_pi", "neg_ka", "neg_pr"}}, "Labeled array for n-sigma values for TPC, TOF, ITS for pions, kaons, protons (positive and negative)"};
  Configurable<bool> cfgUseRun3V2PID{"cfgUseRun3V2PID", true, "True if PID cuts to be used are similar to Run3 v2 PID analysis"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::vector<std::vector<std::shared_ptr<TProfile2D>>> subSample;
  TRandom3* funRndm = new TRandom3(0);

  // Filter command***********
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtLower) && (aod::track::pt < cfgCutPtUpper) && (requireGlobalTrackInFilter()) && (aod::track::tpcChi2NCl < cfgCutTpcChi2NCl) && (aod::track::itsChi2NCl < cfgCutItsChi2NCl) && (nabs(aod::track::dcaZ) < cfgCutTrackDcaZ);

  // Filtering collisions and tracks***********
  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs, aod::Mults>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullEl, aod::pidTOFFullEl>>;

  std::array<float, 6> tofNsigmaCut;
  std::array<float, 6> itsNsigmaCut;
  std::array<float, 6> tpcNsigmaCut;

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
    // Get nSigma values of TPC, TOF, ITS for various particles from the matrix "nSigmas"
    tpcNsigmaCut[kPionUpCut] = nSigmas->getData()[kTPC][kPionUpCut];
    tpcNsigmaCut[kKaonUpCut] = nSigmas->getData()[kTPC][kKaonUpCut];
    tpcNsigmaCut[kProtonUpCut] = nSigmas->getData()[kTPC][kProtonUpCut];
    tpcNsigmaCut[kPionLowCut] = nSigmas->getData()[kTPC][kPionLowCut];
    tpcNsigmaCut[kKaonLowCut] = nSigmas->getData()[kTPC][kKaonLowCut];
    tpcNsigmaCut[kProtonLowCut] = nSigmas->getData()[kTPC][kProtonLowCut];

    tofNsigmaCut[kPionUpCut] = nSigmas->getData()[kTOF][kPionUpCut];
    tofNsigmaCut[kKaonUpCut] = nSigmas->getData()[kTOF][kKaonUpCut];
    tofNsigmaCut[kProtonUpCut] = nSigmas->getData()[kTOF][kProtonUpCut];
    tofNsigmaCut[kPionLowCut] = nSigmas->getData()[kTOF][kPionLowCut];
    tofNsigmaCut[kKaonLowCut] = nSigmas->getData()[kTOF][kKaonLowCut];
    tofNsigmaCut[kProtonLowCut] = nSigmas->getData()[kTOF][kProtonLowCut];

    itsNsigmaCut[kPionUpCut] = nSigmas->getData()[kITS][kPionUpCut];
    itsNsigmaCut[kKaonUpCut] = nSigmas->getData()[kITS][kKaonUpCut];
    itsNsigmaCut[kProtonUpCut] = nSigmas->getData()[kITS][kProtonUpCut];
    itsNsigmaCut[kPionLowCut] = nSigmas->getData()[kITS][kPionLowCut];
    itsNsigmaCut[kKaonLowCut] = nSigmas->getData()[kITS][kKaonLowCut];
    itsNsigmaCut[kProtonLowCut] = nSigmas->getData()[kITS][kProtonLowCut];

    // Define axes
    std::vector<double> ptBin = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
    AxisSpec ptAxis = {ptBin, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec noAxis = {1, 0, 1, "no axis"};
    AxisSpec nSigmaAxis = {200, -5.0, 5.0, "n#sigma"};

    // Add histograms to histogram manager (as in the output object of in AliPhysics)

    // QA hists
    histos.add("hZvtx_after_sel", ";Z (cm)", kTH1F, {{240, -12, 12}});
    histos.add("hCentrality", ";centrality (%)", kTH1F, {{90, 0, 90}});
    histos.add("Hist2D_globalTracks_PVTracks", "", {HistType::kTH2D, {nchAxis, nchAxis}});
    histos.add("Hist2D_cent_nch", "", {HistType::kTH2D, {nchAxis, centAxis}});
    histos.add("hP", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.2, 4.}});
    histos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hPhi", ";#phi", kTH1F, {{100, 0., o2::constants::math::TwoPI}});
    histos.add("hEta", ";#eta", kTH1F, {{100, -2.01, 2.01}});
    histos.add("hDcaXY", ";#it{dca}_{XY}", kTH1F, {{1000, -5, 5}});
    histos.add("hDcaZ", ";#it{dca}_{Z}", kTH1F, {{1000, -5, 5}});
    histos.add("hMeanPt", "", kTProfile, {centAxis});

    // 2D histograms of nSigma
    // before cut
    histos.add("h2DnsigmaPionTpcVsPtBeforeCut", "2D hist of nSigmaTPC vs. pT (pion)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaKaonTpcVsPtBeforeCut", "2D hist of nSigmaTPC vs. pT (kaon)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaProtonTpcVsPtBeforeCut", "2D hist of nSigmaTPC vs. pT (proton)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaPionTofVsPtBeforeCut", "2D hist of nSigmaTOF vs. pT (pion)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaKaonTofVsPtBeforeCut", "2D hist of nSigmaTOF vs. pT (kaon)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaProtonTofVsPtBeforeCut", "2D hist of nSigmaTOF vs. pT (proton)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaPionTpcVsTofBeforeCut", "2D hist of nSigmaTPC vs. nSigmaTOF (pion)", kTH2F, {nSigmaAxis, nSigmaAxis});
    histos.add("h2DnsigmaKaonTpcVsTofBeforeCut", "2D hist of nSigmaTPC vs. nSigmaTOF (kaon)", kTH2F, {nSigmaAxis, nSigmaAxis});
    histos.add("h2DnsigmaProtonTpcVsTofBeforeCut", "2D hist of nSigmaTPC vs. nSigmaTOF (proton)", kTH2F, {nSigmaAxis, nSigmaAxis});
    // after cut
    histos.add("h2DnsigmaPionTpcVsPtAfterCut", "2D hist of nSigmaTPC vs. pT (pion)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaKaonTpcVsPtAfterCut", "2D hist of nSigmaTPC vs. pT (kaon)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaProtonTpcVsPtAfterCut", "2D hist of nSigmaTPC vs. pT (proton)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaPionTofVsPtAfterCut", "2D hist of nSigmaTOF vs. pT (pion)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaKaonTofVsPtAfterCut", "2D hist of nSigmaTOF vs. pT (kaon)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaProtonTofVsPtAfterCut", "2D hist of nSigmaTOF vs. pT (proton)", kTH2F, {ptAxis, nSigmaAxis});
    histos.add("h2DnsigmaPionTpcVsTofAfterCut", "2D hist of nSigmaTPC vs. nSigmaTOF (pion)", kTH2F, {nSigmaAxis, nSigmaAxis});
    histos.add("h2DnsigmaKaonTpcVsTofAfterCut", "2D hist of nSigmaTPC vs. nSigmaTOF (kaon)", kTH2F, {nSigmaAxis, nSigmaAxis});
    histos.add("h2DnsigmaProtonTpcVsTofAfterCut", "2D hist of nSigmaTPC vs. nSigmaTOF (proton)", kTH2F, {nSigmaAxis, nSigmaAxis});

    // Analysis profiles

    histos.add("Prof_A_had", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_C_had", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_D_had", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Bone_had", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Btwo_had", "", {HistType::kTProfile2D, {centAxis, noAxis}});

    histos.add("Prof_A_pi", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_C_pi", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_D_pi", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Bone_pi", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Btwo_pi", "", {HistType::kTProfile2D, {centAxis, noAxis}});

    histos.add("Prof_A_ka", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_C_ka", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_D_ka", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Bone_ka", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Btwo_ka", "", {HistType::kTProfile2D, {centAxis, noAxis}});

    histos.add("Prof_A_prot", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_C_prot", "", {HistType::kTProfile2D, {centAxis, ptAxis}});
    histos.add("Prof_D_prot", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Bone_prot", "", {HistType::kTProfile2D, {centAxis, noAxis}});
    histos.add("Prof_Btwo_prot", "", {HistType::kTProfile2D, {centAxis, noAxis}});

    // initial array
    subSample.resize(cfgNSubsample);
    for (int i = 0; i < cfgNSubsample; i++) {
      subSample[i].resize(20);
    }
    for (int i = 0; i < cfgNSubsample; i++) {
      subSample[i][0] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_A_had", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][1] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_C_had", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][2] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_D_had", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][3] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Bone_had", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][4] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Btwo_had", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));

      subSample[i][5] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_A_pi", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][6] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_C_pi", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][7] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_D_pi", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][8] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Bone_pi", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][9] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Btwo_pi", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));

      subSample[i][10] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_A_ka", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][11] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_C_ka", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][12] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_D_ka", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][13] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Bone_ka", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][14] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Btwo_ka", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));

      subSample[i][15] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_A_prot", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][16] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_C_prot", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      subSample[i][17] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_D_prot", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][18] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Bone_prot", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      subSample[i][19] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("subSample_%d/Prof_Btwo_prot", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
    }
  } // end init

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  template <typename T>
  bool selectionProton(const T& candidate)
  {
    if (!candidate.hasTPC())
      return false;
    int flag = 0; //! pid check main flag

    if (candidate.pt() > cfgCutPtLower && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC) {
        flag = 1;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < cfgnSigmaCutTOF) {
        flag = 1;
      }
    }
    if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < cfgCutPtUpperPID) {
      const float combNSigmaPr = std::sqrt(std::pow(candidate.tpcNSigmaPr(), 2.0) + std::pow(candidate.tofNSigmaPr(), 2.0));
      const float combNSigmaPi = std::sqrt(std::pow(candidate.tpcNSigmaPi(), 2.0) + std::pow(candidate.tofNSigmaPi(), 2.0));
      const float combNSigmaKa = std::sqrt(std::pow(candidate.tpcNSigmaKa(), 2.0) + std::pow(candidate.tofNSigmaKa(), 2.0));

      int flag2 = 0;
      if (combNSigmaPr < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaPi < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaKa < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (!(flag2 > 1) && !(combNSigmaPr > combNSigmaPi) && !(combNSigmaPr > combNSigmaKa)) {
        if (combNSigmaPr < cfgnSigmaCutCombTPCTOF) {
          flag = 1;
        }
      }
    }
    if (flag == 1)
      return true;
    else
      return false;
  }

  template <typename T>
  bool selectionPion(const T& candidate)
  {
    if (!candidate.hasTPC())
      return false;
    int flag = 0; //! pid check main flag

    if (candidate.pt() > cfgCutPtLower && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPC) {
        flag = 1;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaPi()) < cfgnSigmaCutTOF) {
        flag = 1;
      }
    }
    if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < cfgCutPtUpperPID) {
      const float combNSigmaPr = std::sqrt(std::pow(candidate.tpcNSigmaPr(), 2.0) + std::pow(candidate.tofNSigmaPr(), 2.0));
      const float combNSigmaPi = std::sqrt(std::pow(candidate.tpcNSigmaPi(), 2.0) + std::pow(candidate.tofNSigmaPi(), 2.0));
      const float combNSigmaKa = std::sqrt(std::pow(candidate.tpcNSigmaKa(), 2.0) + std::pow(candidate.tofNSigmaKa(), 2.0));

      int flag2 = 0;
      if (combNSigmaPr < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaPi < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaKa < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (!(flag2 > 1) && !(combNSigmaPi > combNSigmaPr) && !(combNSigmaPi > combNSigmaKa)) {
        if (combNSigmaPi < cfgnSigmaCutCombTPCTOF) {
          flag = 1;
        }
      }
    }
    if (flag == 1)
      return true;
    else
      return false;
  }

  template <typename T>
  bool selectionKaon(const T& candidate)
  {
    if (!candidate.hasTPC())
      return false;
    int flag = 0; //! pid check main flag

    if (candidate.pt() > cfgCutPtLower && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPC) {
        flag = 1;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaKa()) < cfgnSigmaCutTOF) {
        flag = 1;
      }
    }
    if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < cfgCutPtUpperPID) {
      const float combNSigmaPr = std::sqrt(std::pow(candidate.tpcNSigmaPr(), 2.0) + std::pow(candidate.tofNSigmaPr(), 2.0));
      const float combNSigmaPi = std::sqrt(std::pow(candidate.tpcNSigmaPi(), 2.0) + std::pow(candidate.tofNSigmaPi(), 2.0));
      const float combNSigmaKa = std::sqrt(std::pow(candidate.tpcNSigmaKa(), 2.0) + std::pow(candidate.tofNSigmaKa(), 2.0));

      int flag2 = 0;
      if (combNSigmaPr < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaPi < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaKa < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (!(flag2 > 1) && !(combNSigmaKa > combNSigmaPi) && !(combNSigmaKa > combNSigmaPr)) {
        if (combNSigmaKa < cfgnSigmaCutCombTPCTOF) {
          flag = 1;
        }
      }
    }
    if (flag == 1)
      return true;
    else
      return false;
  }

  template <typename TTrack>
  int getNsigmaPID(TTrack track)
  {
    // Computing Nsigma arrays for pion, kaon, and protons
    std::array<float, 3> nSigmaTPC = {track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
    std::array<float, 3> nSigmaTOF = {track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr()};
    std::array<float, 3> nSigmaITS = {itsResponse.nSigmaITS<o2::track::PID::Pion>(track), itsResponse.nSigmaITS<o2::track::PID::Kaon>(track), itsResponse.nSigmaITS<o2::track::PID::Proton>(track)};
    int pid = 0; // 0 = not identified, 1 = pion, 2 = kaon, 3 = proton

    std::array<float, 3> nSigmaToUse = cfgUseItsPID ? nSigmaITS : nSigmaTPC;             // Choose which nSigma to use: TPC or ITS
    std::array<float, 6> detectorNsigmaCut = cfgUseItsPID ? itsNsigmaCut : tpcNsigmaCut; // Choose which nSigma to use: TPC or ITS

    bool isPion, isKaon, isProton;
    bool isDetectedPion = nSigmaToUse[PIONS] < detectorNsigmaCut[kPionUpCut] && nSigmaToUse[PIONS] > detectorNsigmaCut[kPionLowCut];
    bool isDetectedKaon = nSigmaToUse[KAONS] < detectorNsigmaCut[kKaonUpCut] && nSigmaToUse[KAONS] > detectorNsigmaCut[kKaonLowCut];
    bool isDetectedProton = nSigmaToUse[PROTONS] < detectorNsigmaCut[kProtonUpCut] && nSigmaToUse[PROTONS] > detectorNsigmaCut[kProtonLowCut];

    bool isTofPion = nSigmaTOF[PIONS] < tofNsigmaCut[kPionUpCut] && nSigmaTOF[PIONS] > tofNsigmaCut[kPionLowCut];
    bool isTofKaon = nSigmaTOF[KAONS] < tofNsigmaCut[kKaonUpCut] && nSigmaTOF[KAONS] > tofNsigmaCut[kKaonLowCut];
    bool isTofProton = nSigmaTOF[PROTONS] < tofNsigmaCut[kProtonUpCut] && nSigmaTOF[PROTONS] > tofNsigmaCut[kProtonLowCut];

    if (track.pt() > cfgPtCutTOF && !track.hasTOF()) {
      return 0;
    } else if (track.pt() > cfgPtCutTOF && track.hasTOF()) {
      isPion = isTofPion && isDetectedPion;
      isKaon = isTofKaon && isDetectedKaon;
      isProton = isTofProton && isDetectedProton;
    } else {
      isPion = isDetectedPion;
      isKaon = isDetectedKaon;
      isProton = isDetectedProton;
    }

    if ((isPion && isKaon) || (isPion && isProton) || (isKaon && isProton)) {
      return 0; // more than one particle satisfy the criteria
    }

    if (isPion) {
      pid = PIONS + 1;
    } else if (isKaon) {
      pid = KAONS + 1;
    } else if (isProton) {
      pid = PROTONS + 1;
    } else {
      return 0; // no particle satisfies the criteria
    }

    return pid; // 0 = not identified, 1 = pion, 2 = kaon, 3 = proton
  }

  // process Data
  void process(AodCollisions::iterator const& coll, aod::BCsWithTimestamps const&, AodTracks const& inputTracks)
  {
    if (!coll.sel8()) {
      return;
    }
    if (cfgUseGoodITSLayerAllCut && !(coll.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))) {
      return;
    }
    if (cfgEvSelkNoSameBunchPileup && !(coll.selection_bit(o2::aod::evsel::kNoSameBunchPileup))) {
      return;
    }
    if (cfgEvSelkNoITSROFrameBorder && !(coll.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    if (cfgEvSelkNoTimeFrameBorder && !(coll.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))) {
      return;
    }
    if (cfgEvSelUseGoodZvtxFT0vsPV && !(coll.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }

    // Centrality
    double cent = 0.0;
    if (cfgCentralityChoice == kFT0C)
      cent = coll.centFT0C();
    else if (cfgCentralityChoice == kFT0A)
      cent = coll.centFT0A();
    else if (cfgCentralityChoice == kFT0M)
      cent = coll.centFT0M();
    else if (cfgCentralityChoice == kFV0A)
      cent = coll.centFV0A();

    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), inputTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), inputTracks.size(), cent);

    // Analysis variables
    int nbinsHad = 20;
    int nbinsPid = 18;
    double binsarray[21] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
    TH1D* fPtProfileHad = new TH1D("fPtProfileHad", "fPtProfileHad", 20, binsarray);
    TH1D* fPtProfilePi = new TH1D("fPtProfilePi", "fPtProfilePi", 20, binsarray);
    TH1D* fPtProfileKa = new TH1D("fPtProfileKa", "fPtProfileKa", 20, binsarray);
    TH1D* fPtProfileProt = new TH1D("fPtProfileProt", "fPtProfileProt", 20, binsarray);
    double pTsumEtaLeftHad = 0.0;
    double nSumEtaLeftHad = 0.0;
    double pTsumEtaRightHad = 0.0;
    double nSumEtaRightHad = 0.0;
    double nSumEtaLeftPi = 0.0;
    double nSumEtaLeftKa = 0.0;
    double nSumEtaLeftProt = 0.0;

    for (const auto& track : inputTracks) { // Loop over tracks

      if (!track.has_collision()) {
        continue;
      }

      if (!track.isPVContributor()) {
        continue;
      }

      if (!(track.itsNCls() > cfgITScluster) || !(track.tpcNClsFound() >= cfgTPCcluster) || !(track.tpcNClsCrossedRows() >= cfgTPCnCrossedRows)) {
        continue;
      }

      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
      histos.fill(HIST("hEta"), track.eta());
      histos.fill(HIST("hPhi"), track.phi());
      histos.fill(HIST("hDcaXY"), track.dcaXY());
      histos.fill(HIST("hDcaZ"), track.dcaZ());

      double trkPt = track.pt();
      double trkEta = track.eta();

      // inclusive charged particles
      if (track.sign() != 0) {
        if (trkEta < cfgCutEtaLeft) {
          fPtProfileHad->Fill(trkPt);
          pTsumEtaLeftHad += trkPt;
          nSumEtaLeftHad += 1.0;
        }
        if (trkEta > cfgCutEtaRight) {
          pTsumEtaRightHad += trkPt;
          nSumEtaRightHad += 1.0;
        }
      }

      // PID QAs before selection
      double nSigmaTpcPi = track.tpcNSigmaPi();
      double nSigmaTpcKa = track.tpcNSigmaKa();
      double nSigmaTpcProt = track.tpcNSigmaPr();
      double nSigmaTofPi = track.tofNSigmaPi();
      double nSigmaTofKa = track.tofNSigmaKa();
      double nSigmaTofProt = track.tofNSigmaPr();
      histos.fill(HIST("h2DnsigmaPionTpcVsPtBeforeCut"), trkPt, nSigmaTpcPi);
      histos.fill(HIST("h2DnsigmaKaonTpcVsPtBeforeCut"), trkPt, nSigmaTpcKa);
      histos.fill(HIST("h2DnsigmaProtonTpcVsPtBeforeCut"), trkPt, nSigmaTpcProt);
      histos.fill(HIST("h2DnsigmaPionTofVsPtBeforeCut"), trkPt, nSigmaTofPi);
      histos.fill(HIST("h2DnsigmaKaonTofVsPtBeforeCut"), trkPt, nSigmaTofKa);
      histos.fill(HIST("h2DnsigmaProtonTofVsPtBeforeCut"), trkPt, nSigmaTofProt);
      histos.fill(HIST("h2DnsigmaPionTpcVsTofBeforeCut"), nSigmaTpcPi, nSigmaTofPi);
      histos.fill(HIST("h2DnsigmaKaonTpcVsTofBeforeCut"), nSigmaTpcKa, nSigmaTofKa);
      histos.fill(HIST("h2DnsigmaProtonTpcVsTofBeforeCut"), nSigmaTpcProt, nSigmaTofProt);

      // identified particles selection
      bool isPion = false;
      bool isKaon = false;
      bool isProton = false;

      if (cfgUseRun3V2PID) {
        int pidVal = getNsigmaPID(track);
        if (pidVal == PIONS + 1)
          isPion = true;
        if (pidVal == KAONS + 1)
          isKaon = true;
        if (pidVal == PROTONS + 1)
          isProton = true;
      } else {
        isPion = selectionPion(track);
        isKaon = selectionKaon(track);
        isProton = selectionProton(track);
      }

      // PID QAs after selection
      if (isPion) {
        histos.fill(HIST("h2DnsigmaPionTpcVsPtAfterCut"), trkPt, nSigmaTpcPi);
        histos.fill(HIST("h2DnsigmaPionTofVsPtAfterCut"), trkPt, nSigmaTofPi);
        histos.fill(HIST("h2DnsigmaPionTpcVsTofAfterCut"), nSigmaTpcPi, nSigmaTofPi);
      }
      if (isKaon) {
        histos.fill(HIST("h2DnsigmaKaonTpcVsPtAfterCut"), trkPt, nSigmaTpcKa);
        histos.fill(HIST("h2DnsigmaKaonTofVsPtAfterCut"), trkPt, nSigmaTofKa);
        histos.fill(HIST("h2DnsigmaKaonTpcVsTofAfterCut"), nSigmaTpcKa, nSigmaTofKa);
      }
      if (isProton) {
        histos.fill(HIST("h2DnsigmaProtonTpcVsPtAfterCut"), trkPt, nSigmaTpcProt);
        histos.fill(HIST("h2DnsigmaProtonTofVsPtAfterCut"), trkPt, nSigmaTofProt);
        histos.fill(HIST("h2DnsigmaProtonTpcVsTofAfterCut"), nSigmaTpcProt, nSigmaTofProt);
      }

      if (track.sign() != 0) {
        if (trkPt < cfgCutPtUpperPID) {
          if (trkEta < cfgCutEtaLeft) {
            if (isPion) {
              fPtProfilePi->Fill(trkPt);
              nSumEtaLeftPi += 1.0;
            }
            if (isKaon) {
              fPtProfileKa->Fill(trkPt);
              nSumEtaLeftKa += 1.0;
            }
            if (isProton && trkPt > cfgCutPtLowerProt) {
              fPtProfileProt->Fill(trkPt);
              nSumEtaLeftProt += 1.0;
            }
          }
        }
      }

    } // End track loop

    // selecting subsample and filling profiles
    float lRandom = funRndm->Rndm();
    int sampleIndex = static_cast<int>(cfgNSubsample * lRandom);

    if (nSumEtaRightHad > 0 && nSumEtaLeftHad > 0) {
      for (int i = 0; i < nbinsHad; i++) {
        histos.get<TProfile2D>(HIST("Prof_A_had"))->Fill(cent, fPtProfileHad->GetBinCenter(i + 1), (fPtProfileHad->GetBinContent(i + 1) / nSumEtaLeftHad));
        histos.get<TProfile2D>(HIST("Prof_C_had"))->Fill(cent, fPtProfileHad->GetBinCenter(i + 1), ((fPtProfileHad->GetBinContent(i + 1) / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));
        histos.get<TProfile2D>(HIST("Prof_Bone_had"))->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        histos.get<TProfile2D>(HIST("Prof_Btwo_had"))->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
        histos.get<TProfile2D>(HIST("Prof_D_had"))->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));

        subSample[sampleIndex][0]->Fill(cent, fPtProfileHad->GetBinCenter(i + 1), (fPtProfileHad->GetBinContent(i + 1) / nSumEtaLeftHad));
        subSample[sampleIndex][1]->Fill(cent, fPtProfileHad->GetBinCenter(i + 1), ((fPtProfileHad->GetBinContent(i + 1) / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][2]->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][3]->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        subSample[sampleIndex][4]->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
      }
    }

    if (nSumEtaRightHad > 0 && nSumEtaLeftHad > 0 && nSumEtaLeftPi > 0) {
      for (int i = 0; i < nbinsPid; i++) {
        histos.get<TProfile2D>(HIST("Prof_A_pi"))->Fill(cent, fPtProfilePi->GetBinCenter(i + 1), (fPtProfilePi->GetBinContent(i + 1) / nSumEtaLeftPi));
        histos.get<TProfile2D>(HIST("Prof_C_pi"))->Fill(cent, fPtProfilePi->GetBinCenter(i + 1), ((fPtProfilePi->GetBinContent(i + 1) / nSumEtaLeftPi) * (pTsumEtaRightHad / nSumEtaRightHad)));
        histos.get<TProfile2D>(HIST("Prof_Bone_pi"))->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        histos.get<TProfile2D>(HIST("Prof_Btwo_pi"))->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
        histos.get<TProfile2D>(HIST("Prof_D_pi"))->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));

        subSample[sampleIndex][5]->Fill(cent, fPtProfilePi->GetBinCenter(i + 1), (fPtProfilePi->GetBinContent(i + 1) / nSumEtaLeftPi));
        subSample[sampleIndex][6]->Fill(cent, fPtProfilePi->GetBinCenter(i + 1), ((fPtProfilePi->GetBinContent(i + 1) / nSumEtaLeftPi) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][7]->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][8]->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        subSample[sampleIndex][9]->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
      }
    }

    if (nSumEtaRightHad > 0 && nSumEtaLeftHad > 0 && nSumEtaLeftKa > 0) {
      for (int i = 0; i < nbinsPid; i++) {
        histos.get<TProfile2D>(HIST("Prof_A_ka"))->Fill(cent, fPtProfileKa->GetBinCenter(i + 1), (fPtProfileKa->GetBinContent(i + 1) / nSumEtaLeftKa));
        histos.get<TProfile2D>(HIST("Prof_C_ka"))->Fill(cent, fPtProfileKa->GetBinCenter(i + 1), ((fPtProfileKa->GetBinContent(i + 1) / nSumEtaLeftKa) * (pTsumEtaRightHad / nSumEtaRightHad)));
        histos.get<TProfile2D>(HIST("Prof_Bone_ka"))->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        histos.get<TProfile2D>(HIST("Prof_Btwo_ka"))->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
        histos.get<TProfile2D>(HIST("Prof_D_ka"))->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));

        subSample[sampleIndex][10]->Fill(cent, fPtProfileKa->GetBinCenter(i + 1), (fPtProfileKa->GetBinContent(i + 1) / nSumEtaLeftKa));
        subSample[sampleIndex][11]->Fill(cent, fPtProfileKa->GetBinCenter(i + 1), ((fPtProfileKa->GetBinContent(i + 1) / nSumEtaLeftKa) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][12]->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][13]->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        subSample[sampleIndex][14]->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
      }
    }
    if (nSumEtaRightHad > 0 && nSumEtaLeftHad > 0 && nSumEtaLeftProt > 0) {
      for (int i = 1; i < nbinsPid; i++) {
        histos.get<TProfile2D>(HIST("Prof_A_prot"))->Fill(cent, fPtProfileProt->GetBinCenter(i + 1), (fPtProfileProt->GetBinContent(i + 1) / nSumEtaLeftProt));
        histos.get<TProfile2D>(HIST("Prof_C_prot"))->Fill(cent, fPtProfileProt->GetBinCenter(i + 1), ((fPtProfileProt->GetBinContent(i + 1) / nSumEtaLeftProt) * (pTsumEtaRightHad / nSumEtaRightHad)));
        histos.get<TProfile2D>(HIST("Prof_Bone_prot"))->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        histos.get<TProfile2D>(HIST("Prof_Btwo_prot"))->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
        histos.get<TProfile2D>(HIST("Prof_D_prot"))->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));

        subSample[sampleIndex][15]->Fill(cent, fPtProfileProt->GetBinCenter(i + 1), (fPtProfileProt->GetBinContent(i + 1) / nSumEtaLeftProt));
        subSample[sampleIndex][16]->Fill(cent, fPtProfileProt->GetBinCenter(i + 1), ((fPtProfileProt->GetBinContent(i + 1) / nSumEtaLeftProt) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][17]->Fill(cent, 0.5, ((pTsumEtaLeftHad / nSumEtaLeftHad) * (pTsumEtaRightHad / nSumEtaRightHad)));
        subSample[sampleIndex][18]->Fill(cent, 0.5, (pTsumEtaLeftHad / nSumEtaLeftHad));
        subSample[sampleIndex][19]->Fill(cent, 0.5, (pTsumEtaRightHad / nSumEtaRightHad));
      }
    }

    fPtProfileHad->Delete();
    fPtProfilePi->Delete();
    fPtProfileKa->Delete();
    fPtProfileProt->Delete();

  } // End process loop
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<V0ptHadPiKaProt>(cfgc)};
  return workflow;
}
