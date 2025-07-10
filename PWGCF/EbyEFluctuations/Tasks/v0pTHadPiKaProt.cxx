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

/// \file v0pTHadPiKaProt.cxx
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
#include <TLorentzVector.h>
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

struct V0pTHadPiKaProt {

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutTpcChi2NCl{"cfgCutTpcChi2NCl", 2.5f, "Maximum TPCchi2NCl"};
  Configurable<float> cfgCutItsChi2NCl{"cfgCutItsChi2NCl", 36.0f, "Maximum ITSchi2NCl"};
  Configurable<float> cfgCutTrackDcaZ{"cfgCutTrackDcaZ", 2.0f, "Maximum DcaZ"};
  Configurable<int> cfgITScluster{"cfgITScluster", 1, "Minimum Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 80, "Minimum Number of TPC cluster"};
  Configurable<int> cfgTPCnCrossedRows{"cfgTPCnCrossedRows", 70, "Minimum Number of TPC crossed-rows"};
  Configurable<float> cfgCutPtUpperTPC{"cfgCutPtUpperTPC", 0.6f, "Upper pT cut for PID using TPC only"};
  Configurable<float> cfgnSigmaCutTPC{"cfgnSigmaCutTPC", 2.0f, "PID nSigma cut for TPC"};
  Configurable<float> cfgnSigmaCutTOF{"cfgnSigmaCutTOF", 2.0f, "PID nSigma cut for TOF"};
  Configurable<float> cfgnSigmaCutCombTPCTOF{"cfgnSigmaCutCombTPCTOF", 2.0f, "PID nSigma combined cut for TPC and TOF"};
  ConfigurableAxis nchAxis{"nchAxis", {5000, 0.5, 5000.5}, ""};
  ConfigurableAxis centAxis{"centAxis", {90, 0., 90.}, "Centrality/Multiplicity percentile bining"};
  Configurable<float> cfgCutPtLower{"cfgCutPtLower", 0.2f, "Lower pT cut"};
  Configurable<float> cfgCutPtUpper{"cfgCutPtUpper", 10.0f, "Higher pT cut"};
  Configurable<float> cfgCutPtUpperPID{"cfgCutPtUpperPID", 6.0f, "Higher pT cut"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "absolute Eta cut"};
  Configurable<float> cfgCutEtaLeft{"cfgCutEtaLeft", 0.8f, "Left end of eta gap"};
  Configurable<float> cfgCutEtaRight{"cfgCutEtaRight", 0.8f, "Right end of eta gap"};
  Configurable<int> cfgNSubsample{"cfgNSubsample", 10, "Number of subsamples"};
  Configurable<int> cfgCentralityChoice{"cfgCentralityChoice", 1, "Which centrality estimator? 0-->FT0M, 1-->FT0C"};
  Configurable<bool> cfgEvSelkNoSameBunchPileup{"cfgEvSelkNoSameBunchPileup", true, "Pileup removal"};
  Configurable<bool> cfgUseGoodITSLayerAllCut{"cfgUseGoodITSLayerAllCut", true, "Remove time interval with dead ITS zone"};

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  std::vector<std::vector<std::shared_ptr<TProfile2D>>> Subsample;
  TRandom3* fRndm = new TRandom3(0);

  // Filter command***********
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < 0.8f) && (aod::track::pt > cfgCutPtLower) && (aod::track::pt < cfgCutPtUpper) && (requireGlobalTrackInFilter()) && (aod::track::tpcChi2NCl < cfgCutTpcChi2NCl) && (aod::track::itsChi2NCl < cfgCutItsChi2NCl) && (nabs(aod::track::dcaZ) < cfgCutTrackDcaZ);

  // Filtering collisions and tracks***********
  using aodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs, aod::Mults>>;
  using aodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullEl, aod::pidTOFFullEl>>;

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
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
    histos.add("hPhi", ";#phi", kTH1F, {{100, 0., 2. * M_PI}});
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
    Subsample.resize(cfgNSubsample);
    for (int i = 0; i < cfgNSubsample; i++) {
      Subsample[i].resize(20);
    }
    for (int i = 0; i < cfgNSubsample; i++) {
      Subsample[i][0] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_A_had", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      Subsample[i][1] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_C_had", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      Subsample[i][2] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_D_had", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      Subsample[i][3] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_Bone_had", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      Subsample[i][4] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_Btwo_had", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));

      Subsample[i][5] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_A_pi", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      Subsample[i][6] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_C_pi", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      Subsample[i][7] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_D_pi", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      Subsample[i][8] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_Bone_pi", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      Subsample[i][9] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_Btwo_pi", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));

      Subsample[i][10] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_A_ka", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      Subsample[i][11] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_C_ka", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      Subsample[i][12] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_D_ka", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      Subsample[i][13] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_Bone_ka", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      Subsample[i][14] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_Btwo_ka", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));

      Subsample[i][15] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_A_prot", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      Subsample[i][16] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_C_prot", i), "", {HistType::kTProfile2D, {centAxis, ptAxis}}));
      Subsample[i][17] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_D_prot", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      Subsample[i][18] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_Bone_prot", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
      Subsample[i][19] = std::get<std::shared_ptr<TProfile2D>>(histos.add(Form("Subsample_%d/Prof_Btwo_prot", i), "", {HistType::kTProfile2D, {centAxis, noAxis}}));
    }
  } // end init

  //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  template <typename T>
  bool selectionProton(const T& candidate)
  {
    if (!candidate.hasTPC())
      return false;
    int flag = 0; //! pid check main flag

    if (candidate.pt() > 0.2f && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC) {
        flag = 1;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < cfgnSigmaCutTOF) {
        flag = 1;
      }
    }
    if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < 5.0f) {
      const float combNSigmaPr = std::sqrt(std::pow(candidate.tpcNSigmaPr(), 2.0) + std::pow(candidate.tofNSigmaPr(), 2.0));
      const float combNSigmaPi = std::sqrt(std::pow(candidate.tpcNSigmaPi(), 2.0) + std::pow(candidate.tofNSigmaPi(), 2.0));
      const float combNSigmaKa = std::sqrt(std::pow(candidate.tpcNSigmaKa(), 2.0) + std::pow(candidate.tofNSigmaKa(), 2.0));

      int flag2 = 0;
      if (combNSigmaPr < 3.0)
        flag2 += 1;
      if (combNSigmaPi < 3.0)
        flag2 += 1;
      if (combNSigmaKa < 3.0)
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

    if (candidate.pt() > 0.2f && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPC) {
        flag = 1;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaPi()) < cfgnSigmaCutTOF) {
        flag = 1;
      }
    }
    if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < 5.0f) {
      const float combNSigmaPr = std::sqrt(std::pow(candidate.tpcNSigmaPr(), 2.0) + std::pow(candidate.tofNSigmaPr(), 2.0));
      const float combNSigmaPi = std::sqrt(std::pow(candidate.tpcNSigmaPi(), 2.0) + std::pow(candidate.tofNSigmaPi(), 2.0));
      const float combNSigmaKa = std::sqrt(std::pow(candidate.tpcNSigmaKa(), 2.0) + std::pow(candidate.tofNSigmaKa(), 2.0));

      int flag2 = 0;
      if (combNSigmaPr < 3.0)
        flag2 += 1;
      if (combNSigmaPi < 3.0)
        flag2 += 1;
      if (combNSigmaKa < 3.0)
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

    if (candidate.pt() > 0.2f && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPC) {
        flag = 1;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaKa()) < cfgnSigmaCutTOF) {
        flag = 1;
      }
    }
    if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < 5.0f) {
      const float combNSigmaPr = std::sqrt(std::pow(candidate.tpcNSigmaPr(), 2.0) + std::pow(candidate.tofNSigmaPr(), 2.0));
      const float combNSigmaPi = std::sqrt(std::pow(candidate.tpcNSigmaPi(), 2.0) + std::pow(candidate.tofNSigmaPi(), 2.0));
      const float combNSigmaKa = std::sqrt(std::pow(candidate.tpcNSigmaKa(), 2.0) + std::pow(candidate.tofNSigmaKa(), 2.0));

      int flag2 = 0;
      if (combNSigmaPr < 3.0)
        flag2 += 1;
      if (combNSigmaPi < 3.0)
        flag2 += 1;
      if (combNSigmaKa < 3.0)
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

  // process Data
  void process(aodCollisions::iterator const& coll, aod::BCsWithTimestamps const&, aodTracks const& inputTracks)
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

    // Centrality
    double cent = 0.0;
    if (cfgCentralityChoice == 0)
      cent = coll.centFT0M();
    else if (cfgCentralityChoice == 1)
      cent = coll.centFT0M();

    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), inputTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), inputTracks.size(), cent);

    // Analysis variables
    double binsarray[21] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
    TH1D* fPt_profile_had = new TH1D("fPt_profile_had", "fPt_profile_had", 20, binsarray);
    TH1D* fPt_profile_pi = new TH1D("fPt_profile_pi", "fPt_profile_pi", 20, binsarray);
    TH1D* fPt_profile_ka = new TH1D("fPt_profile_ka", "fPt_profile_ka", 20, binsarray);
    TH1D* fPt_profile_prot = new TH1D("fPt_profile_prot", "fPt_profile_prot", 20, binsarray);
    double pT_sum_etaLeft_had = 0.0;
    double N_sum_etaLeft_had = 0.0;
    double pT_sum_etaRight_had = 0.0;
    double N_sum_etaRight_had = 0.0;
    double N_sum_etaLeft_pi = 0.0;
    double N_sum_etaLeft_ka = 0.0;
    double N_sum_etaLeft_prot = 0.0;

    for (auto track : inputTracks) { // Loop over tracks

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
          fPt_profile_had->Fill(trkPt);
          pT_sum_etaLeft_had += trkPt;
          N_sum_etaLeft_had += 1.0;
        }
        if (trkEta > cfgCutEtaRight) {
          pT_sum_etaRight_had += trkPt;
          N_sum_etaRight_had += 1.0;
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
      bool isPion = selectionPion(track);
      bool isKaon = selectionKaon(track);
      bool isProton = selectionProton(track);

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
              fPt_profile_pi->Fill(trkPt);
              N_sum_etaLeft_pi += 1.0;
            }
            if (isKaon) {
              fPt_profile_ka->Fill(trkPt);
              N_sum_etaLeft_ka += 1.0;
            }
            if (isProton && trkPt > 0.4) {
              fPt_profile_prot->Fill(trkPt);
              N_sum_etaLeft_prot += 1.0;
            }
          }
        }
      }

    } // End track loop

    // selecting subsample and filling profiles
    float l_Random = fRndm->Rndm();
    int SampleIndex = static_cast<int>(cfgNSubsample * l_Random);

    if (N_sum_etaRight_had > 0 && N_sum_etaLeft_had > 0) {
      for (int i = 0; i < 20; i++) {
        histos.get<TProfile2D>(HIST("Prof_A_had"))->Fill(cent, fPt_profile_had->GetBinCenter(i + 1), (fPt_profile_had->GetBinContent(i + 1) / N_sum_etaLeft_had));
        histos.get<TProfile2D>(HIST("Prof_C_had"))->Fill(cent, fPt_profile_had->GetBinCenter(i + 1), ((fPt_profile_had->GetBinContent(i + 1) / N_sum_etaLeft_had) * (pT_sum_etaRight_had / N_sum_etaRight_had)));
        histos.get<TProfile2D>(HIST("Prof_Bone_had"))->Fill(cent, 0.5, (pT_sum_etaLeft_had / N_sum_etaLeft_had));
        histos.get<TProfile2D>(HIST("Prof_Btwo_had"))->Fill(cent, 0.5, (pT_sum_etaRight_had / N_sum_etaRight_had));
        histos.get<TProfile2D>(HIST("Prof_D_had"))->Fill(cent, 0.5, ((pT_sum_etaLeft_had / N_sum_etaLeft_had) * (pT_sum_etaRight_had / N_sum_etaRight_had)));

        Subsample[SampleIndex][0]->Fill(cent, fPt_profile_had->GetBinCenter(i + 1), (fPt_profile_had->GetBinContent(i + 1) / N_sum_etaLeft_had));
        Subsample[SampleIndex][1]->Fill(cent, fPt_profile_had->GetBinCenter(i + 1), ((fPt_profile_had->GetBinContent(i + 1) / N_sum_etaLeft_had) * (pT_sum_etaRight_had / N_sum_etaRight_had)));
        Subsample[SampleIndex][2]->Fill(cent, 0.5, ((pT_sum_etaLeft_had / N_sum_etaLeft_had) * (pT_sum_etaRight_had / N_sum_etaRight_had)));
        ;
        Subsample[SampleIndex][3]->Fill(cent, 0.5, (pT_sum_etaLeft_had / N_sum_etaLeft_had));
        Subsample[SampleIndex][4]->Fill(cent, 0.5, (pT_sum_etaRight_had / N_sum_etaRight_had));
      }
    }

    if (N_sum_etaRight_had > 0 && N_sum_etaLeft_had > 0 && N_sum_etaLeft_pi > 0) {
      for (int i = 0; i < 18; i++) {
        histos.get<TProfile2D>(HIST("Prof_A_pi"))->Fill(cent, fPt_profile_pi->GetBinCenter(i + 1), (fPt_profile_pi->GetBinContent(i + 1) / N_sum_etaLeft_pi));
        histos.get<TProfile2D>(HIST("Prof_C_pi"))->Fill(cent, fPt_profile_pi->GetBinCenter(i + 1), ((fPt_profile_pi->GetBinContent(i + 1) / N_sum_etaLeft_pi) * (pT_sum_etaRight_had / N_sum_etaRight_had)));
        histos.get<TProfile2D>(HIST("Prof_Bone_pi"))->Fill(cent, 0.5, (pT_sum_etaLeft_had / N_sum_etaLeft_had));
        histos.get<TProfile2D>(HIST("Prof_Btwo_pi"))->Fill(cent, 0.5, (pT_sum_etaRight_had / N_sum_etaRight_had));
        histos.get<TProfile2D>(HIST("Prof_D_pi"))->Fill(cent, 0.5, ((pT_sum_etaLeft_had / N_sum_etaLeft_had) * (pT_sum_etaRight_had / N_sum_etaRight_had)));

        Subsample[SampleIndex][5]->Fill(cent, fPt_profile_pi->GetBinCenter(i + 1), (fPt_profile_pi->GetBinContent(i + 1) / N_sum_etaLeft_pi));
        Subsample[SampleIndex][6]->Fill(cent, fPt_profile_pi->GetBinCenter(i + 1), ((fPt_profile_pi->GetBinContent(i + 1) / N_sum_etaLeft_pi) * (pT_sum_etaRight_had / N_sum_etaRight_had)));
        Subsample[SampleIndex][7]->Fill(cent, 0.5, ((pT_sum_etaLeft_had / N_sum_etaLeft_had) * (pT_sum_etaRight_had / N_sum_etaRight_had)));
        ;
        Subsample[SampleIndex][8]->Fill(cent, 0.5, (pT_sum_etaLeft_had / N_sum_etaLeft_had));
        Subsample[SampleIndex][9]->Fill(cent, 0.5, (pT_sum_etaRight_had / N_sum_etaRight_had));
      }
    }

    if (N_sum_etaRight_had > 0 && N_sum_etaLeft_had > 0 && N_sum_etaLeft_ka > 0) {
      for (int i = 0; i < 18; i++) {
        histos.get<TProfile2D>(HIST("Prof_A_ka"))->Fill(cent, fPt_profile_ka->GetBinCenter(i + 1), (fPt_profile_ka->GetBinContent(i + 1) / N_sum_etaLeft_ka));
        histos.get<TProfile2D>(HIST("Prof_C_ka"))->Fill(cent, fPt_profile_ka->GetBinCenter(i + 1), ((fPt_profile_ka->GetBinContent(i + 1) / N_sum_etaLeft_ka) * (pT_sum_etaRight_had / N_sum_etaRight_had)));
        histos.get<TProfile2D>(HIST("Prof_Bone_ka"))->Fill(cent, 0.5, (pT_sum_etaLeft_had / N_sum_etaLeft_had));
        histos.get<TProfile2D>(HIST("Prof_Btwo_ka"))->Fill(cent, 0.5, (pT_sum_etaRight_had / N_sum_etaRight_had));
        histos.get<TProfile2D>(HIST("Prof_D_ka"))->Fill(cent, 0.5, ((pT_sum_etaLeft_had / N_sum_etaLeft_had) * (pT_sum_etaRight_had / N_sum_etaRight_had)));

        Subsample[SampleIndex][10]->Fill(cent, fPt_profile_ka->GetBinCenter(i + 1), (fPt_profile_ka->GetBinContent(i + 1) / N_sum_etaLeft_ka));
        Subsample[SampleIndex][11]->Fill(cent, fPt_profile_ka->GetBinCenter(i + 1), ((fPt_profile_ka->GetBinContent(i + 1) / N_sum_etaLeft_ka) * (pT_sum_etaRight_had / N_sum_etaRight_had)));
        Subsample[SampleIndex][12]->Fill(cent, 0.5, ((pT_sum_etaLeft_had / N_sum_etaLeft_had) * (pT_sum_etaRight_had / N_sum_etaRight_had)));
        ;
        Subsample[SampleIndex][13]->Fill(cent, 0.5, (pT_sum_etaLeft_had / N_sum_etaLeft_had));
        Subsample[SampleIndex][14]->Fill(cent, 0.5, (pT_sum_etaRight_had / N_sum_etaRight_had));
      }
    }
    if (N_sum_etaRight_had > 0 && N_sum_etaLeft_had > 0 && N_sum_etaLeft_prot > 0) {
      for (int i = 1; i < 18; i++) {
        histos.get<TProfile2D>(HIST("Prof_A_prot"))->Fill(cent, fPt_profile_prot->GetBinCenter(i + 1), (fPt_profile_prot->GetBinContent(i + 1) / N_sum_etaLeft_prot));
        histos.get<TProfile2D>(HIST("Prof_C_prot"))->Fill(cent, fPt_profile_prot->GetBinCenter(i + 1), ((fPt_profile_prot->GetBinContent(i + 1) / N_sum_etaLeft_prot) * (pT_sum_etaRight_had / N_sum_etaRight_had)));
        histos.get<TProfile2D>(HIST("Prof_Bone_prot"))->Fill(cent, 0.5, (pT_sum_etaLeft_had / N_sum_etaLeft_had));
        histos.get<TProfile2D>(HIST("Prof_Btwo_prot"))->Fill(cent, 0.5, (pT_sum_etaRight_had / N_sum_etaRight_had));
        histos.get<TProfile2D>(HIST("Prof_D_prot"))->Fill(cent, 0.5, ((pT_sum_etaLeft_had / N_sum_etaLeft_had) * (pT_sum_etaRight_had / N_sum_etaRight_had)));

        Subsample[SampleIndex][15]->Fill(cent, fPt_profile_prot->GetBinCenter(i + 1), (fPt_profile_prot->GetBinContent(i + 1) / N_sum_etaLeft_prot));
        Subsample[SampleIndex][16]->Fill(cent, fPt_profile_prot->GetBinCenter(i + 1), ((fPt_profile_prot->GetBinContent(i + 1) / N_sum_etaLeft_prot) * (pT_sum_etaRight_had / N_sum_etaRight_had)));
        Subsample[SampleIndex][17]->Fill(cent, 0.5, ((pT_sum_etaLeft_had / N_sum_etaLeft_had) * (pT_sum_etaRight_had / N_sum_etaRight_had)));
        ;
        Subsample[SampleIndex][18]->Fill(cent, 0.5, (pT_sum_etaLeft_had / N_sum_etaLeft_had));
        Subsample[SampleIndex][19]->Fill(cent, 0.5, (pT_sum_etaRight_had / N_sum_etaRight_had));
      }
    }

    fPt_profile_had->Delete();
    fPt_profile_pi->Delete();
    fPt_profile_ka->Delete();
    fPt_profile_prot->Delete();

  } // End process loop
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<V0pTHadPiKaProt>(cfgc)};
  return workflow;
}
