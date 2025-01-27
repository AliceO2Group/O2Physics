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

/// \file netprotonCumulantsMc.cxx
/// \brief Task for analyzing efficiency of proton, and net-proton distributions in MC reconstructed and generated
/// \author Swati Saha

#include <CCDB/BasicCCDBManager.h>
#include <cstdlib>
#include <cmath>
#include <array>
#include <vector>
#include <string>
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/StepTHn.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include <TList.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TRandom3.h>
#include <TMath.h>
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
#include <TF1.h>

namespace o2::aod
{

namespace gen_ebyecolltable
{
DECLARE_SOA_COLUMN(CentralityGen, centralityGen, float);
DECLARE_SOA_COLUMN(NetProtNoGen, netProtNoGen, float);   //! net proton no. in an event
DECLARE_SOA_COLUMN(ProtNoGen, protNoGen, float);         //! proton no. in an event
DECLARE_SOA_COLUMN(AntiProtNoGen, antiProtNoGen, float); //! antiproton no. in an event
} // namespace gen_ebyecolltable

DECLARE_SOA_TABLE(ProtGenCollEbyeTables, "AOD", "PROTGENCOLLEBYETABLE",
                  gen_ebyecolltable::CentralityGen,
                  gen_ebyecolltable::NetProtNoGen,
                  gen_ebyecolltable::ProtNoGen,
                  gen_ebyecolltable::AntiProtNoGen);
using ProtGenCollEbyeTable = ProtGenCollEbyeTables::iterator;

namespace rec_ebyecolltable
{
DECLARE_SOA_COLUMN(CentralityRec, centralityRec, float);
DECLARE_SOA_COLUMN(NetProtNoRec, netProtNoRec, float);   //! net proton no. in an event
DECLARE_SOA_COLUMN(ProtNoRec, protNoRec, float);         //! proton no. in an event
DECLARE_SOA_COLUMN(AntiProtNoRec, antiProtNoRec, float); //! antiproton no. in an event
} // namespace rec_ebyecolltable

DECLARE_SOA_TABLE(ProtRecCollEbyeTables, "AOD", "PROTRECCOLLEBYETABLE",
                  rec_ebyecolltable::CentralityRec,
                  rec_ebyecolltable::NetProtNoRec,
                  rec_ebyecolltable::ProtNoRec,
                  rec_ebyecolltable::AntiProtNoRec);
using ProtRecCollEbyeTable = ProtRecCollEbyeTables::iterator;

DECLARE_SOA_TABLE(ProtRecCollTables, "AOD", "PROTRECCOLLTABLE",
                  o2::soa::Index<>,
                  rec_ebyecolltable::CentralityRec);
using ProtRecCollTable = ProtRecCollTables::iterator;

namespace rec_ebyetracktable
{
DECLARE_SOA_INDEX_COLUMN(ProtRecCollTable, protRecCollTable);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Charge, charge, int);
} // namespace rec_ebyetracktable

DECLARE_SOA_TABLE(ProtRecCompleteEbyeTables, "AOD", "PROTRECCOMPLETEEBYETABLE",
                  o2::soa::Index<>,
                  rec_ebyetracktable::ProtRecCollTableId,
                  rec_ebyetracktable::Pt,
                  rec_ebyetracktable::Eta,
                  rec_ebyetracktable::Charge);
using ProtRecCompleteEbyeTable = ProtRecCompleteEbyeTables::iterator;

} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct NetprotonCumulantsMc {
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // MC
  Configurable<bool> cfgIsMC{"cfgIsMC", true, "Run MC"};
  // tracks
  Configurable<float> cfgCutPtLower{"cfgCutPtLower", 0.2f, "Lower pT cut"};
  Configurable<float> cfgCutPtUpper{"cfgCutPtUpper", 3.0f, "Higher pT cut"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "absolute Eta cut"};
  Configurable<int> cfgPIDchoice{"cfgPIDchoice", 1, "PID selection fucntion choice"};
  Configurable<float> cfgCutPtUpperTPC{"cfgCutPtUpperTPC", 0.6f, "Upper pT cut for PID using TPC only"};
  Configurable<float> cfgnSigmaCutTPC{"cfgnSigmaCutTPC", 2.0f, "PID nSigma cut for TPC"};
  Configurable<float> cfgnSigmaCutTOF{"cfgnSigmaCutTOF", 2.0f, "PID nSigma cut for TOF"};
  Configurable<float> cfgnSigmaCutCombTPCTOF{"cfgnSigmaCutCombTPCTOF", 2.0f, "PID nSigma combined cut for TPC and TOF"};
  Configurable<float> cfgCutTpcChi2NCl{"cfgCutTpcChi2NCl", 2.5f, "Maximum TPCchi2NCl"};
  Configurable<float> cfgCutItsChi2NCl{"cfgCutItsChi2NCl", 36.0f, "Maximum ITSchi2NCl"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};

  // Connect to ccdb
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Filter command for rec (data/MC)***********
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex && (aod::evsel::sel8 == true);
  Filter trackFilter = (nabs(aod::track::eta) < 0.8f) && (aod::track::pt > cfgCutPtLower) && (aod::track::pt < 5.0f) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutTpcChi2NCl) && (aod::track::itsChi2NCl < cfgCutItsChi2NCl) && (aod::track::dcaZ < cfgCutDCAz) && (aod::track::dcaXY < cfgCutDCAxy);

  // filtering collisions and tracks for real data***********
  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullEl, aod::pidTOFFullEl>>;

  // filtering collisions and tracks for MC rec data***********
  using MyMCRecCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFDDMs, aod::McCollisionLabels>>;
  using MyMCRecCollision = MyMCRecCollisions::iterator;
  using MyMCTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullEl, aod::pidTOFFullEl, aod::McTrackLabels>>;
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFV0As, aod::CentFDDMs>;

  // // Filter command for gen ***************
  // Filter mcCollFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;
  // Filter mcParticleFilter = ((aod::mcparticle::pt > cfgCutPtLower) && (aod::mcparticle::pt < 5.0f) && (nabs(aod::mcparticle::eta) < 0.8f));
  // // filtering collisions and particles for MC gen data***********
  // using myMCGenCollision = soa::Filtered<aod::McCollisions>::iterator;
  // using myMCParticles = soa::Filtered<aod::McParticles>;

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
    // Define your axes
    // Constant bin width axis
    AxisSpec vtxZAxis = {100, -20, 20};
    // Variable bin width axis
    std::vector<double> ptBinning = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    std::vector<double> centBining = {0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90};
    AxisSpec centAxis = {centBining, "Multiplicity percentile from FT0M (%)"};
    AxisSpec netprotonAxis = {41, -20.5, 20.5, "net-proton number"};
    AxisSpec protonAxis = {21, -0.5, 20.5, "proton number"};
    AxisSpec antiprotonAxis = {21, -0.5, 20.5, "antiproton number"};

    // histograms for events
    histos.add("hZvtx_after_sel", "Vertex dist. after event selection;Z (cm)", kTH1F, {vtxZAxis});
    histos.add("hCentrec", "MCRec Multiplicity percentile from FT0M (%)", kTH1F, {{100, 0.0, 100.0}});
    // tracks Rec level histograms
    histos.add("hrecPtAll", "Reconstructed All particles;#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hrecPtProton", "Reconstructed Protons;#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hrecPtAntiproton", "Reconstructed Antiprotons;#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hrecPhiAll", "Reconstructed All particles;#phi", kTH1F, {{100, 0., 7.}});
    histos.add("hrecPhiProton", "Reconstructed Protons;#phi", kTH1F, {{100, 0., 7.}});
    histos.add("hrecPhiAntiproton", "Reconstructed Antiprotons;#phi", kTH1F, {{100, 0., 7.}});
    histos.add("hrecEtaAll", "Reconstructed All particles;#eta", kTH1F, {{100, -2.01, 2.01}});
    histos.add("hrecEtaProton", "Reconstructed Proton;#eta", kTH1F, {{100, -2.01, 2.01}});
    histos.add("hrecEtaAntiproton", "Reconstructed Antiprotons;#eta", kTH1F, {{100, -2.01, 2.01}});
    histos.add("hrecPtDistProtonVsCentrality", "Reconstructed proton number vs centrality in 2D", kTH2F, {ptAxis, centAxis});
    histos.add("hrecPtDistAntiprotonVsCentrality", "Reconstructed antiproton number vs centrality in 2D", kTH2F, {ptAxis, centAxis});
    histos.add("hrecNetProtonVsCentrality", "Reconstructed net-proton number vs centrality in 2D", kTH2F, {netprotonAxis, centAxis});
    histos.add("hrecProtonVsCentrality", "Reconstructed proton number vs centrality in 2D", kTH2F, {protonAxis, centAxis});
    histos.add("hrecAntiprotonVsCentrality", "Reconstructed antiproton number vs centrality in 2D", kTH2F, {antiprotonAxis, centAxis});
    histos.add("hrecProfileTotalProton", "Reconstructed total proton number vs. centrality", kTProfile, {centAxis});
    histos.add("hrecProfileProton", "Reconstructed proton number vs. centrality", kTProfile, {centAxis});
    histos.add("hrecProfileAntiproton", "Reconstructed antiproton number vs. centrality", kTProfile, {centAxis});

    if (cfgIsMC) {
      // MC event counts
      histos.add("hMC", "MC Event statistics", kTH1F, {{10, 0.0f, 10.0f}});
      histos.add("hCentgen", "MCGen Multiplicity percentile from FT0M (%)", kTH1F, {{100, 0.0, 100.0}});
      // tracks Gen level histograms
      histos.add("hgenPtAll", "Generated All particles;#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
      histos.add("hgenPtProton", "Generated Protons;#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
      histos.add("hgenPtAntiproton", "Generated Antiprotons;#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
      histos.add("hgenPhiAll", "Generated All particles;#phi", kTH1F, {{100, 0., 7.}});
      histos.add("hgenPhiProton", "Generated Protons;#phi", kTH1F, {{100, 0., 7.}});
      histos.add("hgenPhiAntiproton", "Generated Antiprotons;#phi", kTH1F, {{100, 0., 7.}});
      histos.add("hgenEtaAll", "Generated All particles;#eta", kTH1F, {{100, -2.01, 2.01}});
      histos.add("hgenEtaProton", "Generated Proton;#eta", kTH1F, {{100, -2.01, 2.01}});
      histos.add("hgenEtaAntiproton", "Generated Antiprotons;#eta", kTH1F, {{100, -2.01, 2.01}});
      histos.add("hgenPtDistProtonVsCentrality", "Generated proton number vs centrality in 2D", kTH2F, {ptAxis, centAxis});
      histos.add("hgenPtDistAntiprotonVsCentrality", "Generated antiproton number vs centrality in 2D", kTH2F, {ptAxis, centAxis});
      histos.add("hrecTruePtProton", "Reconstructed pdgcode verified protons;#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
      histos.add("hrecTruePtAntiproton", "Reconstructed pdgcode verified Antiprotons;#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
      histos.add("hgenNetProtonVsCentrality", "Generated net-proton number vs centrality in 2D", kTH2F, {netprotonAxis, centAxis});
      histos.add("hgenProtonVsCentrality", "Generated proton number vs centrality in 2D", kTH2F, {protonAxis, centAxis});
      histos.add("hgenAntiprotonVsCentrality", "Generated antiproton number vs centrality in 2D", kTH2F, {antiprotonAxis, centAxis});
      histos.add("hgenProfileTotalProton", "Generated total proton number vs. centrality", kTProfile, {centAxis});
      histos.add("hgenProfileProton", "Generated proton number vs. centrality", kTProfile, {centAxis});
      histos.add("hgenProfileAntiproton", "Generated antiproton number vs. centrality", kTProfile, {centAxis});
    }
  }

  template <typename T>
  bool selectionPIDold(const T& candidate)
  {
    //! PID checking as done in Run2 my analysis
    //! ----------------------------------------------------------------------
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
  bool selectionPIDnew(const T& candidate)
  {
    //! if pt < threshold
    if (candidate.pt() > 0.2f && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC && std::abs(candidate.tpcNSigmaPi()) > cfgnSigmaCutTPC && std::abs(candidate.tpcNSigmaKa()) > cfgnSigmaCutTPC) {
        return true;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC && std::abs(candidate.tpcNSigmaPi()) > cfgnSigmaCutTPC && std::abs(candidate.tpcNSigmaKa()) > cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < cfgnSigmaCutTOF && std::abs(candidate.tofNSigmaPi()) > cfgnSigmaCutTOF && std::abs(candidate.tofNSigmaKa()) > cfgnSigmaCutTOF) {
        return true;
      }
    }

    //! if pt > threshold
    if (candidate.pt() > cfgCutPtUpperTPC) {
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC && std::abs(candidate.tpcNSigmaPi()) > cfgnSigmaCutTPC && std::abs(candidate.tpcNSigmaKa()) > cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < cfgnSigmaCutTOF && std::abs(candidate.tofNSigmaPi()) > cfgnSigmaCutTOF && std::abs(candidate.tofNSigmaKa()) > cfgnSigmaCutTOF) {
        return true;
      }
    }
    return false;
  }

  Produces<aod::ProtGenCollEbyeTables> genEbyeCollisions; //! MC Gen table creation

  void processMCGen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {
    histos.fill(HIST("hMC"), 0.5);
    if (std::abs(mcCollision.posZ()) < cfgCutVertex) {
      histos.fill(HIST("hMC"), 1.5);
    }
    auto cent = 0;

    int nchInel = 0;
    for (const auto& mcParticle : mcParticles) {
      auto pdgcode = std::abs(mcParticle.pdgCode());
      if (mcParticle.isPhysicalPrimary() && (pdgcode == 211 || pdgcode == 321 || pdgcode == 2212 || pdgcode == 11 || pdgcode == 13)) {
        if (std::abs(mcParticle.eta()) < 1.0) {
          nchInel = nchInel + 1;
        }
      }
    }
    if (nchInel > 0 && std::abs(mcCollision.posZ()) < cfgCutVertex)
      histos.fill(HIST("hMC"), 2.5);
    std::vector<int64_t> selectedEvents(collisions.size());
    int nevts = 0;

    for (const auto& collision : collisions) {
      if (!collision.sel8() || std::abs(collision.mcCollision().posZ()) > cfgCutVertex) {
        continue;
      }
      cent = collision.centFT0M();

      selectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    selectedEvents.resize(nevts);
    const auto evtReconstructedAndSelected = std::find(selectedEvents.begin(), selectedEvents.end(), mcCollision.globalIndex()) != selectedEvents.end();
    histos.fill(HIST("hMC"), 3.5);
    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    histos.fill(HIST("hMC"), 4.5);
    histos.fill(HIST("hCentgen"), cent);

    // creating phi, pt, eta dstribution of generted MC particles

    float nProt = 0.0;
    float nAntiprot = 0.0;

    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary()) {
        if ((mcParticle.pt() > cfgCutPtLower) && (mcParticle.pt() < 5.0f) && (std::abs(mcParticle.eta()) < 0.8f)) {
          histos.fill(HIST("hgenPtAll"), mcParticle.pt());
          histos.fill(HIST("hgenEtaAll"), mcParticle.eta());
          histos.fill(HIST("hgenPhiAll"), mcParticle.phi());

          if (std::abs(mcParticle.pdgCode()) == 2212 /*&& std::abs(mcParticle.y()) < 0.5*/) {
            if (mcParticle.pdgCode() == 2212) {
              histos.fill(HIST("hgenPtProton"), mcParticle.pt()); //! hist for p gen
              histos.fill(HIST("hgenPtDistProtonVsCentrality"), mcParticle.pt(), cent);
              histos.fill(HIST("hgenEtaProton"), mcParticle.eta());
              histos.fill(HIST("hgenPhiProton"), mcParticle.phi());
              if (mcParticle.pt() < cfgCutPtUpper)
                nProt = nProt + 1.0;
            }
            if (mcParticle.pdgCode() == -2212) {
              histos.fill(HIST("hgenPtAntiproton"), mcParticle.pt()); //! hist for anti-p gen
              histos.fill(HIST("hgenPtDistAntiprotonVsCentrality"), mcParticle.pt(), cent);
              histos.fill(HIST("hgenEtaAntiproton"), mcParticle.eta());
              histos.fill(HIST("hgenPhiAntiproton"), mcParticle.phi());
              if (mcParticle.pt() < cfgCutPtUpper)
                nAntiprot = nAntiprot + 1.0;
            }
          }
        }
      }
    } //! end particle loop

    float netProt = nProt - nAntiprot;
    histos.fill(HIST("hgenNetProtonVsCentrality"), netProt, cent);
    histos.fill(HIST("hgenProtonVsCentrality"), nProt, cent);
    histos.fill(HIST("hgenAntiprotonVsCentrality"), nAntiprot, cent);
    histos.fill(HIST("hgenProfileTotalProton"), cent, (nProt + nAntiprot));
    histos.fill(HIST("hgenProfileProton"), cent, nProt);
    histos.fill(HIST("hgenProfileAntiproton"), cent, nAntiprot);
    genEbyeCollisions(cent, netProt, nProt, nAntiprot);
  }
  PROCESS_SWITCH(NetprotonCumulantsMc, processMCGen, "Process Generated", true);

  Produces<aod::ProtRecCollEbyeTables> recEbyeCollisions;             //! MC Rec table creation
  Produces<aod::ProtRecCollTables> recCollisions;                     //! MC Rec table creation
  Produces<aod::ProtRecCompleteEbyeTables> recEbyeCompleteCollisions; //! MC Rec table creation with tracks

  void processMCRec(MyMCRecCollision const& collision, MyMCTracks const& tracks, aod::McCollisions const&, aod::McParticles const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    auto cent = collision.centFT0M();
    histos.fill(HIST("hCentrec"), cent);
    histos.fill(HIST("hMC"), 5.5);
    histos.fill(HIST("hZvtx_after_sel"), collision.posZ());
    recCollisions(cent);

    float nProt = 0.0;
    float nAntiprot = 0.0;

    // Start of the Monte-Carlo reconstructed tracks
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) //! check if track has corresponding MC particle
      {
        continue;
      }
      if (!track.isPVContributor()) //! track check as used in data
      {
        continue;
      }

      auto particle = track.mcParticle();
      if (particle.isPhysicalPrimary()) {
        histos.fill(HIST("hrecPtAll"), particle.pt());
        histos.fill(HIST("hrecEtaAll"), particle.eta());
        histos.fill(HIST("hrecPhiAll"), particle.phi());

        bool trackSelected = false;
        if (cfgPIDchoice == 0)
          trackSelected = selectionPIDold(track);
        if (cfgPIDchoice == 1)
          trackSelected = selectionPIDnew(track);

        if (trackSelected) {
          recEbyeCompleteCollisions(recCollisions.lastIndex(), particle.pt(), particle.eta(), track.sign());
          if (track.sign() > 0) {
            histos.fill(HIST("hrecPtProton"), particle.pt()); //! hist for p rec
            histos.fill(HIST("hrecPtDistProtonVsCentrality"), particle.pt(), cent);
            histos.fill(HIST("hrecEtaProton"), particle.eta());
            histos.fill(HIST("hrecPhiProton"), particle.phi());
            if (particle.pt() < cfgCutPtUpper)
              nProt = nProt + 1.0;
            if (particle.pdgCode() == 2212) {
              histos.fill(HIST("hrecTruePtProton"), particle.pt()); //! hist for p purity
            }
          }
          if (track.sign() < 0) {
            histos.fill(HIST("hrecPtAntiproton"), particle.pt()); //! hist for anti-p rec
            histos.fill(HIST("hrecPtDistAntiprotonVsCentrality"), particle.pt(), cent);
            histos.fill(HIST("hrecEtaAntiproton"), particle.eta());
            histos.fill(HIST("hrecPhiAntiproton"), particle.phi());
            if (particle.pt() < cfgCutPtUpper)
              nAntiprot = nAntiprot + 1.0;
            if (particle.pdgCode() == -2212) {
              histos.fill(HIST("hrecTruePtAntiproton"), particle.pt()); //! hist for anti-p purity
            }
          }
        } //! checking PID
      } //! checking if primary
    } //! end track loop

    float netProt = nProt - nAntiprot;
    histos.fill(HIST("hrecNetProtonVsCentrality"), netProt, cent);
    histos.fill(HIST("hrecProtonVsCentrality"), nProt, cent);
    histos.fill(HIST("hrecAntiprotonVsCentrality"), nAntiprot, cent);
    histos.fill(HIST("hrecProfileTotalProton"), cent, (nProt + nAntiprot));
    histos.fill(HIST("hrecProfileProton"), cent, nProt);
    histos.fill(HIST("hrecProfileAntiproton"), cent, nAntiprot);
    recEbyeCollisions(cent, netProt, nProt, nAntiprot);
  }
  PROCESS_SWITCH(NetprotonCumulantsMc, processMCRec, "Process Generated", true);

  void processDataRec(AodCollisions::iterator const& coll, aod::BCsWithTimestamps const&, AodTracks const& inputTracks)
  {
    if (!coll.sel8()) {
      return;
    }

    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());
    // variables
    auto cent = coll.centFT0M();
    histos.fill(HIST("hCentrec"), cent);
    recCollisions(cent);

    float nProt = 0.0;
    float nAntiprot = 0.0;

    // Start of the Monte-Carlo reconstructed tracks
    for (const auto& track : inputTracks) {
      if (!track.isPVContributor()) //! track check as used in data
      {
        continue;
      }

      histos.fill(HIST("hrecPtAll"), track.pt());
      histos.fill(HIST("hrecEtaAll"), track.eta());
      histos.fill(HIST("hrecPhiAll"), track.phi());

      bool trackSelected = false;
      if (cfgPIDchoice == 0)
        trackSelected = selectionPIDold(track);
      if (cfgPIDchoice == 1)
        trackSelected = selectionPIDnew(track);

      if (trackSelected) {
        recEbyeCompleteCollisions(recCollisions.lastIndex(), track.pt(), track.eta(), track.sign());
        if (track.sign() > 0) {
          histos.fill(HIST("hrecPtProton"), track.pt()); //! hist for p rec
          histos.fill(HIST("hrecPtDistProtonVsCentrality"), track.pt(), cent);
          histos.fill(HIST("hrecEtaProton"), track.eta());
          histos.fill(HIST("hrecPhiProton"), track.phi());
          if (track.pt() < cfgCutPtUpper)
            nProt = nProt + 1.0;
        }
        if (track.sign() < 0) {
          histos.fill(HIST("hrecPtAntiproton"), track.pt()); //! hist for anti-p rec
          histos.fill(HIST("hrecPtDistAntiprotonVsCentrality"), track.pt(), cent);
          histos.fill(HIST("hrecEtaAntiproton"), track.eta());
          histos.fill(HIST("hrecPhiAntiproton"), track.phi());
          if (track.pt() < cfgCutPtUpper)
            nAntiprot = nAntiprot + 1.0;
        }
      } //! checking PID
    } //! end track loop

    float netProt = nProt - nAntiprot;
    histos.fill(HIST("hrecNetProtonVsCentrality"), netProt, cent);
    histos.fill(HIST("hrecProtonVsCentrality"), nProt, cent);
    histos.fill(HIST("hrecAntiprotonVsCentrality"), nAntiprot, cent);
    histos.fill(HIST("hrecProfileTotalProton"), cent, (nProt + nAntiprot));
    histos.fill(HIST("hrecProfileProton"), cent, nProt);
    histos.fill(HIST("hrecProfileAntiproton"), cent, nAntiprot);
    recEbyeCollisions(cent, netProt, nProt, nAntiprot);
  }
  PROCESS_SWITCH(NetprotonCumulantsMc, processDataRec, "Process real data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<NetprotonCumulantsMc>(cfgc)};
  return workflow;
}
