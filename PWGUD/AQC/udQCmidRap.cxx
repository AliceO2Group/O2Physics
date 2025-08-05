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
/// \brief A task for Asynchronus Quality Control for Ultra-perimpheral and Diffraction (AQC-UD) two tracks (pion, kaon, muon and electron) candidates in midrapidity
/// \author Anisa Khatun
/// \author Paul Buehler
/// \since 17.01.2023

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/BCRange.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/DataModel/FT0Corrected.h"
#include "PWGUD/Core/UDHelpers.h"
#include "Framework/StaticFor.h"
#include "TLorentzVector.h"
#include "TMath.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UDQCmid {

  SliceCache cache;
  Preslice<aod::Zdcs> perBCzdc = aod::zdc::bcId;
  Preslice<aod::Calos> perBCcalo = aod::calo::bcId;

  // global variables
  float maxdEdxTPC;
  float maxdEdxTOF;

  // get a DGCutparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  Configurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};
  // Configurable<bool> withAmbTrackAnalysis{"ambiguousTracks", false, "with ambiguous tracks analysis"};
  // Configurable<bool> withAmbFwdTrackAnalysis{"ambiguousFwdTracks", false, "with ambiguous forward tracks analysis"};
  //  Configurable<bool> doCleanFITBC{"doCleanFITBC", false, "Require cleanFIT in compatible BCs"};

  // structures to hold information about the possible BCs the ambiguous tracks/FwdTracks belong to
  o2::dataformats::bcRanges abcrs = o2::dataformats::bcRanges("ambiguous_tracks");
  o2::dataformats::bcRanges afbcrs = o2::dataformats::bcRanges("ambiguous_fwdtracks");

  // inivinitialize HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {}};

  // define abbreviations
  using CCs = soa::Join<aod::Collisions, aod::EvSels>;
  using CC = CCs::iterator;
  using BCs = soa::Join<aod::BCs, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using TCs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::TrackSelectionExtension, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullMu, aod::pidTPCFullEl, aod::TOFSignal, aod::pidTOFbeta>;
  using FWs = aod::FwdTracks;
  using ATs = aod::AmbiguousTracks;
  using AFTs = aod::AmbiguousFwdTracks;

  Partition<TCs> goodTracks = requireGlobalTrackInFilter();

  void init(InitContext& context)
  {
    // initialize global variables
    maxdEdxTPC = 0.;
    maxdEdxTOF = 0.;
    diffCuts = (DGCutparHolder)DGCuts;

    // add histograms for the different process functions
    if (context.mOptions.get<bool>("processMain")) {
      // collisions
      registry.add("collisions/Stat", "Cut statistics; Selection criterion; Collisions", {HistType::kTH1F, {{20, -0.5, 19.5}}});
      registry.add("collisions/Tracks", "Number of tracks; Number of tracks; Collisions", {HistType::kTH1F, {{300, 0.5, 300.5}}});
      registry.add("collisions/vtxTracks", "Number of vertex tracks; Number of contributors; Collisions", {HistType::kTH1F, {{300, 0.5, 300.5}}});
      registry.add("collisions/globalTracks", "Number of global tracks; Number of global tracks; Collisions", {HistType::kTH1F, {{300, 0.5, 300.5}}});
      registry.add("collisions/tResvsrTOFTracks", "Number of PV tracks with TOF hit versus collision time resolution; Collision time resolution [ns]; Fraction of PV tracks with TOF hit; Collisions", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}});
      registry.add("collisions/tResvsTOFTrkNoPV", "Number of No PV tracks with TOF hit versus collision time resolution; Collision time resolution [ns]; Fraction of No PV tracks with TOF hit; Collisions", {HistType::kTH2F, {{1000, 0., 1.E3}, {101, -0.01, 1.01}}});

      // tracks
      registry.add("tracks/dEdxTPC", "TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("tracks/dEdxTOF", "TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});

      // DG
      registry.add("DG/etapt", "DG: eta versus pT of all tracks; eta of track; p_T of track [GeV/c^2]; Tracks", {HistType::kTH2F, {{80, -2., 2.}, {100, 0., 5.}}});
      registry.add("DG/dEdxTPC", "DG: TPC signal versus signed track momentum; Signed track momentum [GeV/c]; TPC signal; Tracks", {HistType::kTH2F, {{120, -6., 6.}, {1000, 0., 1000.}}});
      registry.add("DG/dEdxTOF", "DG: TOF signal versus signed track momentum; Signed track momentum [GeV/c]; TOF signal; Tracks", {HistType::kTH2F, {{1000, 0., 10.}, {1000, 0., 10.}}});

      registry.add("DG/hMassAll", "DG: Invariant mass of pions; Invarian mass [GeV/c^2]", {HistType::kTH1F, {{1000, 0., 10.}}});
      registry.add("DG/hMassFIT", "DG: Invariant mass of pions; Invarian mass [GeV/c^2]", {HistType::kTH1F, {{1000, 0., 10.}}});

      registry.add("DG/hMassZDC", "DG: Invariant mass of pions; Invarian mass [GeV/c^2]", {HistType::kTH1F, {{1000, 0., 10.}}});
      registry.add("DG/hMassTOF", "DG: Invariant mass of pions; Invarian mass [GeV/c^2]", {HistType::kTH1F, {{1000, 0., 10.}}});
      registry.add("DG/hMassFWD", "DG: Invariant mass of pions; Invarian mass [GeV/c^2]", {HistType::kTH1F, {{1000, 0., 10.}}});
      registry.add("DG/hMassGlobalTrk", "DG: Invariant mass of pions; Invarian mass [GeV/c^2]", {HistType::kTH1F, {{1000, 0., 10.}}});
      registry.add("DG/hMassITSTrk", "DG: Invariant mass of pions; Invarian mass [GeV/c^2]", {HistType::kTH1F, {{1000, 0., 10.}}});

      registry.add("DG/hMassAmbigous", "DG: Invariant mass of pions; Invarian mass [GeV/c^2]", {HistType::kTH1F, {{1000, 0., 10.}}});
      registry.add("DG/hMassAmbigousFWD", "DG: Invariant mass of pions; Invarian mass [GeV/c^2]", {HistType::kTH1F, {{1000, 0., 10.}}});

      registry.add("DG/etaphi", "DG: Eta versus Phi; eta ; #phi ", {HistType::kTH2F, {{80, -2., 2.}, {120, 0., 6.28}}});
      registry.add("DG/etaphi1", "DG: Eta versus Phi; eta ; #phi ", {HistType::kTH2F, {{80, -2., 2.}, {120, 0., 6.28}}});
      registry.add("DG/etaphi2", "DG: Eta versus Phi; eta ; #phi ", {HistType::kTH2F, {{80, -2., 2.}, {120, 0., 6.28}}});
      registry.add("DG/etaphi3", "DG: Eta versus Phi; eta ; #phi ", {HistType::kTH2F, {{80, -2., 2.}, {120, 0., 6.28}}});
      registry.add("DG/etaphi4", "DG: Eta versus Phi; eta ; #phi ", {HistType::kTH2F, {{80, -2., 2.}, {120, 0., 6.28}}});
      registry.add("DG/etaphi5", "DG: Eta versus Phi; eta ; #phi ", {HistType::kTH2F, {{80, -2., 2.}, {120, 0., 6.28}}});
      registry.add("DG/etaphi6", "DG: Eta versus Phi; eta ; #phi ", {HistType::kTH2F, {{80, -2., 2.}, {120, 0., 6.28}}});
      registry.add("DG/etaphi7", "DG: Eta versus Phi; eta ; #phi ", {HistType::kTH2F, {{80, -2., 2.}, {120, 0., 6.28}}});
      registry.add("DG/etaphi8", "DG: Eta versus Phi; eta ; #phi ", {HistType::kTH2F, {{80, -2., 2.}, {120, 0., 6.28}}});

      registry.add("DG/IVMptSys2PVtrk", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions 2 PV tracks", {HistType::kTH2F, {{2000, 0., 20.}, {2000, 0., 20.0}}});
      registry.add("DG/IVMptSys2PVtrk1", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions 2 PV tracks", {HistType::kTH2F, {{2000, 0., 20.}, {2000, 0., 20.0}}});
      registry.add("DG/IVMptSys2PVtrk2", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions 2 PV tracks", {HistType::kTH2F, {{2000, 0., 20.}, {2000, 0., 20.0}}});
      registry.add("DG/IVMptSys2PVtrk3", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions 2 PV tracks", {HistType::kTH2F, {{2000, 0., 20.}, {2000, 0., 20.0}}});
      registry.add("DG/IVMptSys2PVtrk4", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions 2 PV tracks", {HistType::kTH2F, {{2000, 0., 20.}, {2000, 0., 20.0}}});
      registry.add("DG/IVMptSys2PVtrk5", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions 2 PV tracks", {HistType::kTH2F, {{2000, 0., 20.}, {2000, 0., 20.0}}});
      registry.add("DG/IVMptSys2PVtrk6", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions 2 PV tracks", {HistType::kTH2F, {{2000, 0., 20.}, {2000, 0., 20.0}}});
      registry.add("DG/IVMptSys2PVtrk7", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions 2 PV tracks", {HistType::kTH2F, {{2000, 0., 20.}, {2000, 0., 20.0}}});
      registry.add("DG/IVMptSys2PVtrk8", "DG: Invariant mass versus p_{T, system}; Invarian mass [GeV/c^2]; p_{T, system} [GeV/c]; DG collisions 2 PV tracks", {HistType::kTH2F, {{2000, 0., 20.}, {2000, 0., 20.0}}});
    }

    if (context.mOptions.get<bool>("processFewProng")) {
      registry.add("fpStat", "#fpStat", {HistType::kTH1F, {{2, 0.5, 2.5}}});
      registry.add("allPVC", "#allPVC", {HistType::kTH1F, {{200, 0.5, 200.5}}});
      registry.add("fpPVC", "#fpPVC", {HistType::kTH1F, {{200, 0.5, 200.5}}});
      registry.add("fpPVC1", "#fpPVC1", {HistType::kTH1F, {{200, 0.5, 200.5}}});
      registry.add("fpPVC2", "#fpPVC2", {HistType::kTH1F, {{200, 0.5, 200.5}}});
    }
  }

  //...............................................................................................................................................
  void processMain(CC const& collision, BCs const& bct0s,
                   TCs const& tracks, FWs const& fwdtracks, ATs const& /*ambtracks*/, AFTs const& /*ambfwdtracks*/,
                   aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/, aod::FDDs const& /*fdds*/,
                   aod::Zdcs& zdcs)
  {
    LOGF(debug, "<UDQCmid. Collision %d", collision.globalIndex());
    LOGF(debug, "<UDQCmid> Start %i", abcrs.size());

    bool isDGcandidate = true;
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(0., isDGcandidate * 1.);

    // update collision histograms
    // tracks
    registry.get<TH1>(HIST("collisions/Tracks"))->Fill(tracks.size());
    // vertex tracks normally gives PV contributors from collisions
    registry.get<TH1>(HIST("collisions/vtxTracks"))->Fill(collision.numContrib());
    // global tracks

    goodTracks.bindTable(tracks);
    registry.get<TH1>(HIST("collisions/globalTracks"))->Fill(goodTracks.size());

    // 12. net charge and invariant mass
    // bool goodetas = false;
    // bool goodpts = false;
    bool ispipiCand = false;
    auto netCharge = 0;
    auto lvtmp = TLorentzVector();
    auto ivm = TLorentzVector();
    if (isDGcandidate) {

      // which particle hypothesis? // Pion hypothesis has been used currently
      auto mass2Use = constants::physics::MassPionCharged;
      if (diffCuts.pidHypothesis() == 321) {
        mass2Use = constants::physics::MassKaonCharged;
      }

      if (diffCuts.pidHypothesis() == 13) {
        mass2Use = constants::physics::MassMuon;
      }

      if (diffCuts.pidHypothesis() == 11) {
        mass2Use = constants::physics::MassElectron;
      }

      // check also pt and eta of tracks
      for (auto const& track : tracks) {

        // define Lorentz vector to create invariant mass
        lvtmp.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), mass2Use);
        LOGF(debug, "mass %f track pt %f/%f eta %f/%f", mass2Use, track.pt(), lvtmp.Perp(), track.eta(), lvtmp.Eta());
        if (track.pt() <= diffCuts.minPt() || track.pt() >= diffCuts.maxPt()) {
          continue;
        }
        if (track.eta() <= diffCuts.minEta() || track.eta() >= diffCuts.maxEta()) {
          continue;
        }
        netCharge += track.sign();
        ivm += lvtmp;
      }

      if (collision.numContrib() == 2) {
        ispipiCand = true;
        for (auto const& track : tracks) {
          if (track.isPVContributor()) {
            if (diffCuts.pidHypothesis() == 211) {
              if (std::abs(track.tpcNSigmaPi()) > diffCuts.maxNSigmaTPC()) {
                ispipiCand = false;
              }
            } else if (diffCuts.pidHypothesis() == 321) {
              if (std::abs(track.tpcNSigmaKa()) > diffCuts.maxNSigmaTPC()) {
                ispipiCand = false;
              }
            } else if (diffCuts.pidHypothesis() == 13) {
              if (std::abs(track.tpcNSigmaMu()) > diffCuts.maxNSigmaTPC()) {
                ispipiCand = false;
              }
            } else if (diffCuts.pidHypothesis() == 11) {
              if (std::abs(track.tpcNSigmaEl()) > diffCuts.maxNSigmaTPC()) {
                ispipiCand = false;
              }
            }
          }
        } // trks
        if (ispipiCand) {
          registry.get<TH2>(HIST("DG/IVMptSys2PVtrk"))->Fill(ivm.M(), ivm.Pt());
          registry.get<TH2>(HIST("DG/etaphi"))->Fill(ivm.Eta(), ivm.Phi());
          if (netCharge == 0)
            registry.get<TH1>(HIST("DG/hMassAll"))->Fill(ivm.M());
        }
      } // coll
    } // dgcand

    // loop over all tracks
    float rgtrwTOF = 0.;
    float norgtrwTOF = 0.;
    for (auto const& track : tracks) {
      // update PV track stats
      if (track.isPVContributor()) {
        // update dEdx histograms
        registry.get<TH2>(HIST("tracks/dEdxTPC"))->Fill(track.tpcInnerParam() / track.sign(), track.tpcSignal());

        if (track.tpcSignal() > maxdEdxTPC) {
          maxdEdxTPC = track.tpcSignal();
          // LOGF(debug, "<UDQCmid> New maxdEdx TPC %f", maxdEdxTPC);
        }

        // TOF hit?
        if (track.hasTOF()) {
          registry.get<TH2>(HIST("tracks/dEdxTOF"))->Fill(track.p() / track.sign(), track.beta());
          if (track.tofSignal() > maxdEdxTOF) {
            maxdEdxTOF = track.tofSignal();
            // LOGF(debug, "<UDQCmid> New maxdEdx TOF %f", maxdEdxTOF);
          }

          // No vertex track with TOF hit?
          if (!track.isPVContributor()) {
            norgtrwTOF += 1.;
          }

          // vertex track with TOF hit?
          if (track.isPVContributor()) {
            rgtrwTOF += 1.;
          }
        }
      }
    } // closing track loop

    // fraction of No PV tracks with TOF hit
    if (collision.numContrib() > 0) {
      norgtrwTOF /= collision.numContrib();
    }

    // fraction of PV tracks with TOF hit
    if (collision.numContrib() > 0) {
      rgtrwTOF /= collision.numContrib();
    }
    // LOGF(debug, "<UDQCmid> PV tracks with TOF: %f [1]", rgtrwTOF);
    registry.get<TH2>(HIST("collisions/tResvsrTOFTracks"))->Fill(collision.collisionTimeRes(), rgtrwTOF);
    registry.get<TH2>(HIST("collisions/tResvsTOFTrkNoPV"))->Fill(collision.collisionTimeRes(), norgtrwTOF);

    // is it a DG candidate?
    // 1. DG = no FIT signal in compatible BCs
    // 2. & no ZDC signal in compatible BCs
    // 3. & number of forward tracks = 0
    // 4. Check for global tracks which are no vtx tracks
    // 5. check a given bc for possible ambiguous Tracks
    // 6. check a given bc for possible ambiguous FwdTracks
    // 7. fraction of PV tracks with TOF hit
    isDGcandidate = true;

    // get BCrange to test for FIT signals
    auto bcSlice = udhelpers::compatibleBCs(collision, diffCuts.NDtcoll(), bct0s, diffCuts.minNBCs());

    // 1. no FIT signal in bcSlice / collision
    for (auto const& bc : bcSlice) {
      if (udhelpers::FITveto(bc, diffCuts)) {
        isDGcandidate = false;
        break;
      }
    }
    /* if (doCleanFITBC) {
       for (auto const& bc : bcSlice) {
         if (!udhelpers::cleanFIT(bc, diffCuts.maxFITtime(), diffCuts.FITAmpLimits())) {
           isDGcandidate = false;
           break;
         }
       }
     }*/
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(1., isDGcandidate * 1.);

    // Invariant mass with 2 PV contributors and all contributors
    if (isDGcandidate) {
      if (ispipiCand) {
        registry.get<TH2>(HIST("DG/IVMptSys2PVtrk1"))->Fill(ivm.M(), ivm.Pt());
        registry.get<TH2>(HIST("DG/etaphi1"))->Fill(ivm.Eta(), ivm.Phi());
        if ((ivm.Pt() < 0.2) && (netCharge == 0)) {
          registry.get<TH1>(HIST("DG/hMassFIT"))->Fill(ivm.M());
        }
      }

      for (auto const& track : tracks) {
        if (track.isPVContributor()) {
          registry.get<TH2>(HIST("DG/etapt"))->Fill(track.eta(), track.pt(), 1.);
          LOGF(debug, "dEdx TPC %f TOF %i %f", track.tpcSignal(), track.hasTOF(), track.hasTOF() ? track.tofSignal() : 0.);
          registry.get<TH2>(HIST("DG/dEdxTPC"))->Fill(track.tpcInnerParam() / track.sign(), track.tpcSignal());

          if (track.hasTOF()) {
            registry.get<TH2>(HIST("DG/dEdxTOF"))->Fill(track.p() / track.sign(), track.beta());
          } // fill TOF
        } // pv contributor
      }
    } // Inavariant mass after FIT

    // 2. no Zdc signal in bcSlice
    bool isZDCcandidate = isDGcandidate;
    std::vector<float> lims(10, 0.);
    for (auto const& bc : bcSlice) {
      if (!udhelpers::cleanZDC(bc, zdcs, lims, cache)) {
        isZDCcandidate = false;
        break;
      }
    }
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(2., isZDCcandidate * 1.);

    if (isZDCcandidate) {
      if (ispipiCand) {
        registry.get<TH2>(HIST("DG/IVMptSys2PVtrk2"))->Fill(ivm.M(), ivm.Pt());
        registry.get<TH2>(HIST("DG/etaphi2"))->Fill(ivm.Eta(), ivm.Phi());
        if ((ivm.Pt() < 0.2) && (netCharge == 0)) {
          registry.get<TH1>(HIST("DG/hMassZDC"))->Fill(ivm.M());
        }
      }
    } // Invariant mass after ZDC

    // 3. number of forward tracks = 0
    isDGcandidate &= (fwdtracks.size() == 0);
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(3., isDGcandidate * 1.);

    if (isDGcandidate) {
      if (ispipiCand) {
        registry.get<TH2>(HIST("DG/IVMptSys2PVtrk3"))->Fill(ivm.M(), ivm.Pt());
        registry.get<TH2>(HIST("DG/etaphi3"))->Fill(ivm.Eta(), ivm.Phi());
        if ((ivm.Pt() < 0.2) && (netCharge == 0)) {
          registry.get<TH1>(HIST("DG/hMassFWD"))->Fill(ivm.M());
        }
      }
    }

    // 4. Check for global tracks
    bool globalAndVtx = isDGcandidate;
    for (auto const& track : tracks) {
      if (track.isGlobalTrack() && track.isPVContributor()) {
        globalAndVtx = false;
      }
    }
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(4., globalAndVtx * 1.);

    if (globalAndVtx) {
      if (ispipiCand) {
        registry.get<TH2>(HIST("DG/IVMptSys2PVtrk4"))->Fill(ivm.M(), ivm.Pt());
        registry.get<TH2>(HIST("DG/etaphi4"))->Fill(ivm.Eta(), ivm.Phi());
        if ((ivm.Pt() < 0.2) && (netCharge == 0)) {
          registry.get<TH1>(HIST("DG/hMassGlobalTrk"))->Fill(ivm.M());
        }
      }
    }

    // 5. Check for ITS only tracks
    bool ITSOnlytrks = isDGcandidate;
    for (auto const& track : tracks) {
      if (track.hasITS() && track.isPVContributor() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF()) {
        ITSOnlytrks = false;
      }
    }
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(5., ITSOnlytrks * 1.);

    if (ITSOnlytrks) {
      if (ispipiCand) {
        registry.get<TH2>(HIST("DG/IVMptSys2PVtrk5"))->Fill(ivm.M(), ivm.Pt());
        registry.get<TH2>(HIST("DG/etaphi5"))->Fill(ivm.Eta(), ivm.Phi());
        if ((ivm.Pt() < 0.2) && (netCharge == 0)) {
          registry.get<TH1>(HIST("DG/hMassITSTrk"))->Fill(ivm.M());
        }
      }
    }

    // 6. check a given bc for possible ambiguous Tracks
    auto noAmbTracks = isDGcandidate;
    for (auto const& bc : bcSlice) {
      if (abcrs.isInRange(bc.globalIndex())) {
        noAmbTracks = false;
        break;
      }
    }

    registry.get<TH1>(HIST("collisions/Stat"))->Fill(6., noAmbTracks * 1.); // noAmbTracks

    if (noAmbTracks) {
      if (ispipiCand) {
        registry.get<TH2>(HIST("DG/IVMptSys2PVtrk6"))->Fill(ivm.M(), ivm.Pt());
        registry.get<TH2>(HIST("DG/etaphi6"))->Fill(ivm.Eta(), ivm.Phi());
        if ((ivm.Pt() < 0.2) && (netCharge == 0)) {
          registry.get<TH1>(HIST("DG/hMassAmbigous"))->Fill(ivm.M());
        }
      }
    }

    // 7. check a given bc for possible ambiguous FwdTracks
    auto noAmbFwdTracks = isDGcandidate;
    for (auto const& bc : bcSlice) {
      if (afbcrs.isInRange(bc.globalIndex())) {
        noAmbFwdTracks = false;
        break;
      }
    }

    registry.get<TH1>(HIST("collisions/Stat"))->Fill(7., noAmbFwdTracks * 1.); // noAmbFwdTracks

    if (noAmbFwdTracks) {
      if (ispipiCand) {
        registry.get<TH2>(HIST("DG/IVMptSys2PVtrk7"))->Fill(ivm.M(), ivm.Pt());
        registry.get<TH2>(HIST("DG/etaphi7"))->Fill(ivm.Eta(), ivm.Phi());
        if ((ivm.Pt() < 0.2) && (netCharge == 0)) {
          registry.get<TH1>(HIST("DG/hMassAmbigousFWD"))->Fill(ivm.M());
        }
      }
    }

    // 8. fraction of PV tracks with TOF hit
    isDGcandidate &= (rgtrwTOF >= diffCuts.minRgtrwTOF());
    registry.get<TH1>(HIST("collisions/Stat"))->Fill(8., isDGcandidate * 1.);

    if (isDGcandidate) {
      if (ispipiCand) {
        registry.get<TH2>(HIST("DG/IVMptSys2PVtrk8"))->Fill(ivm.M(), ivm.Pt());
        registry.get<TH2>(HIST("DG/etaphi8"))->Fill(ivm.Eta(), ivm.Phi());
        if ((ivm.Pt() < 0.2) && (netCharge == 0)) {
          registry.get<TH1>(HIST("DG/hMassTOF"))->Fill(ivm.M());
        }
      }
    }
  }
  PROCESS_SWITCH(UDQCmid, processMain, "Process Main", true);

  //.....................................................................................................................
  // Distribution of number of PV contributors for all collisions and those with empty FT0
  void processFewProng(CC const& collision, BCs const& /*bct0s*/,
                       aod::FT0s const& /*ft0s*/, aod::FV0As const& /*fv0as*/, aod::FDDs const& /*fdds*/)
  {
    // count collisions
    registry.get<TH1>(HIST("fpStat"))->Fill(1., 1.);
    registry.get<TH1>(HIST("allPVC"))->Fill(collision.numContrib(), 1.);

    // check FT0 to be empty
    auto bc = collision.foundBC_as<BCs>();
    if (udhelpers::cleanFT0(bc, diffCuts.maxFITtime(), 0., 0.)) {
      // only collisions with empty FT0 arrive here
      registry.get<TH1>(HIST("fpStat"))->Fill(2., 1.);

      // update #PV contributors in collisions with empty FT0
      registry.get<TH1>(HIST("fpPVC"))->Fill(collision.numContrib(), 1.);

      if (udhelpers::cleanFV0(bc, diffCuts.maxFITtime(), 0.)) {
        // only collisions with empty FV0 and FT0 arrive here
        registry.get<TH1>(HIST("fpStat"))->Fill(3., 1.);

        // update #PV contributors in collisions with empty FT0 && FV0
        registry.get<TH1>(HIST("fpPVC1"))->Fill(collision.numContrib(), 1.);

        if (udhelpers::cleanFDD(bc, diffCuts.maxFITtime(), 0., 0.)) {
          // only collisions with empty FV0 arrive here
          registry.get<TH1>(HIST("fpStat"))->Fill(4., 1.);

          // update #PV contributors in collisions with empty FT0 && FV0&& FDCC
          registry.get<TH1>(HIST("fpPVC2"))->Fill(collision.numContrib(), 1.);
        } // fdd
      } // fvo

    } // ft0
  }
  PROCESS_SWITCH(UDQCmid, processFewProng, "Process FewProng", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UDQCmid>(cfgc, TaskName{"udQCmidRap"}),
  };
}
