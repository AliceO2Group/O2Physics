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

#include <algorithm>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "TPDGCode.h"
#include "RecoDecay.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct ThreePartCorr {
  Service<o2::ccdb::BasicCCDBManager> CCDB;

  // Histogram registry
  HistogramRegistry MECorrRegistry{"MECorrRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry SECorrRegistry{"SECorrRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry MCRegistry{"MCRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry QARegistry{"QARegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  // Collision & Event filters
  Filter CollCent = aod::cent::centFT0C > 0.0f && aod::cent::centFT0C < 90.0f;
  Filter CollZvtx = nabs(aod::collision::posZ) < 7.0f;
  Filter MCCollZvtx = nabs(aod::mccollision::posZ) < 7.0f;
  Filter EvSelect = aod::evsel::sel8 == true;

  // V0 filters
  Filter V0Pt = aod::v0data::pt > 0.6f && aod::v0data::pt < 12.0f;
  Filter V0Eta = nabs(aod::v0data::eta) < 0.72f;

  // Track filters
  Filter TrackPt = aod::track::pt > 0.2f && aod::track::pt < 3.0f;
  Filter TrackEta = nabs(aod::track::eta) < 0.8f;
  Filter GlobalTracks = requireGlobalTrackInFilter();

  // Particle filters
  Filter ParticlePt = aod::mcparticle::pt > 0.2f && aod::mcparticle::pt < 3.0f;
  Filter ParticleEta = nabs(aod::mcparticle::eta) < 0.8f;

  // Table aliases - Data
  using MyFilteredCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>>;
  using MyFilteredCollision = MyFilteredCollisions::iterator;
  using MyFilteredV0s = soa::Filtered<aod::V0Datas>;
  using MyFilteredTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
                                                   aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr,
                                                   aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFbeta>>;

  // Table aliases - MC
  using MyFilteredMCGenCollision = soa::Filtered<aod::McCollisions>::iterator;
  using MyFilteredMCParticles = soa::Filtered<aod::McParticles>;
  using MyFilteredMCRecCollision = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>>::iterator;
  using MyFilteredMCTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::McTrackLabels>>;

  // Mixed-events binning policy
  SliceCache cache;
  ConfigurableAxis ConfCentBins{"ConfCentBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f}, "ME Centrality binning"};
  ConfigurableAxis ConfZvtxBins{"ConfZvtxBins", {VARIABLE_WIDTH, -7.0f, -5.0f, -3.0f, -1.0f, 0.0f, 1.0f, 3.0f, 5.0f, 7.0f}, "ME Zvtx binning"};
  using BinningType = ColumnBinningPolicy<aod::cent::CentFT0C, aod::collision::PosZ>;

  BinningType CollBinning{{ConfCentBins, ConfZvtxBins}, true};
  Pair<MyFilteredCollisions, MyFilteredV0s, MyFilteredTracks, BinningType> pair{CollBinning, 5, -1, &cache};

  // Process configurables
  Configurable<bool> ConfFilterSwitch{"ConfFilterSwitch", false, "Switch for the FakeV0Filter function"};

  // Particle masses
  Double_t massLambda = o2::constants::physics::MassLambda0;
  Double_t DGaussSigma = 0.0021;

  // Efficiency histograms
  TH1D** hEffPions = new TH1D*[2];
  TH1D** hEffKaons = new TH1D*[2];
  TH1D** hEffProtons = new TH1D*[2];

  // Correlation variables
  Int_t T_Sign;
  Double_t CandMass;
  Double_t* A_PID;

  Double_t DeltaPhi, DeltaEta;

  //================================================================================================================================================================================================================

  void init(InitContext const&)
  {

    const AxisSpec CentralityAxis{ConfCentBins};
    const AxisSpec ZvtxAxis{ConfZvtxBins};
    const AxisSpec PhiAxis{36, (-1. / 2) * M_PI, (3. / 2) * M_PI};
    const AxisSpec EtaAxis{32, -1.52, 1.52};
    const AxisSpec V0PtAxis{114, 0.6, 12};
    const AxisSpec TrackPtAxis{28, 0.2, 3};
    const AxisSpec LambdaInvMassAxis{100, 1.08, 1.16};

    QARegistry.add("hTrackPt", "hTrackPt", {HistType::kTH1D, {{100, 0, 4}}});
    QARegistry.add("hTrackEta", "hTrackEta", {HistType::kTH1D, {{100, -1, 1}}});
    QARegistry.add("hTrackPhi", "hTrackPhi", {HistType::kTH1D, {{100, (-1. / 2) * M_PI, (5. / 2) * M_PI}}});
    QARegistry.add("hEventCentrality", "hEventCentrality", {HistType::kTH1D, {{CentralityAxis}}});
    QARegistry.add("hEventZvtx", "hEventZvtx", {HistType::kTH1D, {{ZvtxAxis}}});

    QARegistry.add("hdEdx", "hdEdx", {HistType::kTH2D, {{56, 0.2, 3.0}, {180, 20, 200}}});
    QARegistry.add("hdEdxPion", "hdEdxPion", {HistType::kTH2D, {{56, 0.2, 3.0}, {180, 20, 200}}});
    QARegistry.add("hdEdxKaon", "hdEdxKaon", {HistType::kTH2D, {{56, 0.2, 3.0}, {180, 20, 200}}});
    QARegistry.add("hdEdxProton", "hdEdxProton", {HistType::kTH2D, {{56, 0.2, 3.0}, {180, 20, 200}}});
    QARegistry.add("hBeta", "hBeta", {HistType::kTH2D, {{56, 0.2, 3.0}, {70, 0.4, 1.1}}});
    QARegistry.add("hBetaPion", "hBetaPion", {HistType::kTH2D, {{56, 0.2, 3.0}, {70, 0.4, 1.1}}});
    QARegistry.add("hBetaKaon", "hBetaKaon", {HistType::kTH2D, {{56, 0.2, 3.0}, {70, 0.4, 1.1}}});
    QARegistry.add("hBetaProton", "hBetaProton", {HistType::kTH2D, {{56, 0.2, 3.0}, {70, 0.4, 1.1}}});
    QARegistry.add("hNSigmaPion", "hNSigmaPion", {HistType::kTH2D, {{201, -5.025, 5.025}, {201, -5.025, 5.025}}});
    QARegistry.add("hNSigmaKaon", "hNSigmaKaon", {HistType::kTH2D, {{201, -5.025, 5.025}, {201, -5.025, 5.025}}});
    QARegistry.add("hNSigmaProton", "hNSigmaProton", {HistType::kTH2D, {{201, -5.025, 5.025}, {201, -5.025, 5.025}}});

    QARegistry.add("hInvMassLambda", "hInvMassLambda", {HistType::kTH3D, {{LambdaInvMassAxis}, {V0PtAxis}, {CentralityAxis}}});
    QARegistry.add("hInvMassAntiLambda", "hInvMassAntiLambda", {HistType::kTH3D, {{LambdaInvMassAxis}, {V0PtAxis}, {CentralityAxis}}});

    MCRegistry.add("hGenPionP", "hGenPionP", {HistType::kTH1D, {TrackPtAxis}});
    MCRegistry.add("hGenPionN", "hGenPionN", {HistType::kTH1D, {TrackPtAxis}});
    MCRegistry.add("hGenKaonP", "hGenKaonP", {HistType::kTH1D, {TrackPtAxis}});
    MCRegistry.add("hGenKaonN", "hGenKaonN", {HistType::kTH1D, {TrackPtAxis}});
    MCRegistry.add("hGenProtonP", "hGenProtonP", {HistType::kTH1D, {TrackPtAxis}});
    MCRegistry.add("hGenProtonN", "hGenProtonN", {HistType::kTH1D, {TrackPtAxis}});
    MCRegistry.add("hRecPionP", "hRecPionP", {HistType::kTH1D, {TrackPtAxis}});
    MCRegistry.add("hRecPionN", "hRecPionN", {HistType::kTH1D, {TrackPtAxis}});
    MCRegistry.add("hRecKaonP", "hRecKaonP", {HistType::kTH1D, {TrackPtAxis}});
    MCRegistry.add("hRecKaonN", "hRecKaonN", {HistType::kTH1D, {TrackPtAxis}});
    MCRegistry.add("hRecProtonP", "hRecProtonP", {HistType::kTH1D, {TrackPtAxis}});
    MCRegistry.add("hRecProtonN", "hRecProtonN", {HistType::kTH1D, {TrackPtAxis}});

    SECorrRegistry.add("hSameLambdaPion_SGNL", "Same-event #Lambda - #pi correlator (SGNL region)", {HistType::kTHnSparseD, {{PhiAxis}, {EtaAxis}, {CentralityAxis}, {ZvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    SECorrRegistry.add("hSameLambdaPion_SB", "Same-event #Lambda - #pi correlator (SB region)", {HistType::kTHnSparseD, {{PhiAxis}, {EtaAxis}, {CentralityAxis}, {ZvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    SECorrRegistry.add("hSameLambdaKaon_SGNL", "Same-event #Lambda - K correlator (SGNL region)", {HistType::kTHnSparseD, {{PhiAxis}, {EtaAxis}, {CentralityAxis}, {ZvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    SECorrRegistry.add("hSameLambdaKaon_SB", "Same-event #Lambda - K correlator (SB region)", {HistType::kTHnSparseD, {{PhiAxis}, {EtaAxis}, {CentralityAxis}, {ZvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    SECorrRegistry.add("hSameLambdaProton_SGNL", "Same-event #Lambda - p correlator (SGNL region)", {HistType::kTHnSparseD, {{PhiAxis}, {EtaAxis}, {CentralityAxis}, {ZvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    SECorrRegistry.add("hSameLambdaProton_SB", "Same-event #Lambda - p correlator (SB region)", {HistType::kTHnSparseD, {{PhiAxis}, {EtaAxis}, {CentralityAxis}, {ZvtxAxis}, {2, -2, 2}, {2, -2, 2}}});

    MECorrRegistry.add("hMixLambdaPion_SGNL", "Mixed-event #Lambda - #pi correlator (SGNL region)", {HistType::kTHnSparseD, {{PhiAxis}, {EtaAxis}, {CentralityAxis}, {ZvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    MECorrRegistry.add("hMixLambdaPion_SB", "Mixed-event #Lambda - #pi correlator (SB region)", {HistType::kTHnSparseD, {{PhiAxis}, {EtaAxis}, {CentralityAxis}, {ZvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    MECorrRegistry.add("hMixLambdaKaon_SGNL", "Mixed-event #Lambda - K correlator (SGNL region)", {HistType::kTHnSparseD, {{PhiAxis}, {EtaAxis}, {CentralityAxis}, {ZvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    MECorrRegistry.add("hMixLambdaKaon_SB", "Mixed-event #Lambda - K correlator (SB region)", {HistType::kTHnSparseD, {{PhiAxis}, {EtaAxis}, {CentralityAxis}, {ZvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    MECorrRegistry.add("hMixLambdaProton_SGNL", "Mixed-event #Lambda - p correlator (SGNL region)", {HistType::kTHnSparseD, {{PhiAxis}, {EtaAxis}, {CentralityAxis}, {ZvtxAxis}, {2, -2, 2}, {2, -2, 2}}});
    MECorrRegistry.add("hMixLambdaProton_SB", "Mixed-event #Lambda - p correlator (SB region)", {HistType::kTHnSparseD, {{PhiAxis}, {EtaAxis}, {CentralityAxis}, {ZvtxAxis}, {2, -2, 2}, {2, -2, 2}}});

    CCDB->setURL("http://alice-ccdb.cern.ch");
    CCDB->setCaching(true);
    TList* EfficiencyList = CCDB->getForTimeStamp<TList>("Users/j/jstaa/Efficiency/ChargedParticles", 1);
    hEffPions[0] = static_cast<TH1D*>(EfficiencyList->FindObject("hEfficiencyPionP"));
    hEffPions[1] = static_cast<TH1D*>(EfficiencyList->FindObject("hEfficiencyPionN"));
    hEffKaons[0] = static_cast<TH1D*>(EfficiencyList->FindObject("hEfficiencyKaonP"));
    hEffKaons[1] = static_cast<TH1D*>(EfficiencyList->FindObject("hEfficiencyKaonN"));
    hEffProtons[0] = static_cast<TH1D*>(EfficiencyList->FindObject("hEfficiencyProtonP"));
    hEffProtons[1] = static_cast<TH1D*>(EfficiencyList->FindObject("hEfficiencyProtonN"));
  }

  //================================================================================================================================================================================================================

  void processSame(MyFilteredCollision const& collision, MyFilteredV0s const& v0s, MyFilteredTracks const& tracks)
  {

    QARegistry.fill(HIST("hEventCentrality"), collision.centFT0C());
    QARegistry.fill(HIST("hEventZvtx"), collision.posZ());

    // Start of the Track QA
    for (const auto& track : tracks) {
      A_PID = TrackPID(track);
      if (A_PID[1] < 4.0) {
        QARegistry.fill(HIST("hTrackPt"), track.pt());
        QARegistry.fill(HIST("hTrackEta"), track.eta());
        QARegistry.fill(HIST("hTrackPhi"), track.phi());
        QARegistry.fill(HIST("hdEdx"), track.p(), track.tpcSignal());
        QARegistry.fill(HIST("hBeta"), track.p(), track.beta());
        if (A_PID[0] == 0.0) { // Pions
          QARegistry.fill(HIST("hdEdxPion"), track.p(), track.tpcSignal());
          QARegistry.fill(HIST("hBetaPion"), track.p(), track.beta());
          QARegistry.fill(HIST("hNSigmaPion"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        } else if (A_PID[0] == 1.0) { // Kaons
          QARegistry.fill(HIST("hdEdxKaon"), track.p(), track.tpcSignal());
          QARegistry.fill(HIST("hBetaKaon"), track.p(), track.beta());
          QARegistry.fill(HIST("hNSigmaKaon"), track.tpcNSigmaKa(), track.tofNSigmaKa());
        } else if (A_PID[0] == 2.0) { // Protons
          QARegistry.fill(HIST("hdEdxProton"), track.p(), track.tpcSignal());
          QARegistry.fill(HIST("hBetaProton"), track.p(), track.beta());
          QARegistry.fill(HIST("hNSigmaProton"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        }
      }
    }
    // End of the Track QA

    // Start of the V0-Track Correlations
    for (const auto& trigger : v0s) {
      if (V0Filters(trigger)) {

        T_Sign = V0Sign(trigger);
        if (T_Sign == 1) {
          CandMass = trigger.mLambda();
          QARegistry.fill(HIST("hInvMassLambda"), trigger.mLambda(), trigger.pt(), collision.centFT0C());
        } else if (T_Sign == -1) {
          CandMass = trigger.mAntiLambda();
          QARegistry.fill(HIST("hInvMassAntiLambda"), trigger.mAntiLambda(), trigger.pt(), collision.centFT0C());
        }

        for (const auto& associate : tracks) {
          if (TrackFilters(trigger, associate)) {
            if (FakeV0Filter(trigger, associate)) {

              A_PID = TrackPID(associate);
              DeltaPhi = DeltaPhiShift(trigger.phi(), associate.phi());
              DeltaEta = trigger.eta() - associate.eta();

              if (CandMass >= massLambda - 4 * DGaussSigma && CandMass <= massLambda + 4 * DGaussSigma) {
                if (A_PID[0] == 0.0) { // Pions
                  SECorrRegistry.fill(HIST("hSameLambdaPion_SGNL"), DeltaPhi, DeltaEta, collision.centFT0C(), collision.posZ(), T_Sign, associate.sign(), 1. / TrackEff(hEffPions, associate.sign(), associate.pt()));
                } else if (A_PID[0] == 1.0) { // Kaons
                  SECorrRegistry.fill(HIST("hSameLambdaKaon_SGNL"), DeltaPhi, DeltaEta, collision.centFT0C(), collision.posZ(), T_Sign, associate.sign(), 1. / TrackEff(hEffKaons, associate.sign(), associate.pt()));
                } else if (A_PID[0] == 2.0) { // Protons
                  SECorrRegistry.fill(HIST("hSameLambdaProton_SGNL"), DeltaPhi, DeltaEta, collision.centFT0C(), collision.posZ(), T_Sign, associate.sign(), 1. / TrackEff(hEffProtons, associate.sign(), associate.pt()));
                }
              } else if (CandMass >= massLambda - 8 * DGaussSigma && CandMass <= massLambda + 8 * DGaussSigma) {
                if (A_PID[0] == 0.0) { // Pions
                  SECorrRegistry.fill(HIST("hSameLambdaPion_SB"), DeltaPhi, DeltaEta, collision.centFT0C(), collision.posZ(), T_Sign, associate.sign(), 1. / TrackEff(hEffPions, associate.sign(), associate.pt()));
                } else if (A_PID[0] == 1.0) { // Kaons
                  SECorrRegistry.fill(HIST("hSameLambdaKaon_SB"), DeltaPhi, DeltaEta, collision.centFT0C(), collision.posZ(), T_Sign, associate.sign(), 1. / TrackEff(hEffKaons, associate.sign(), associate.pt()));
                } else if (A_PID[0] == 2.0) { // Protons
                  SECorrRegistry.fill(HIST("hSameLambdaProton_SB"), DeltaPhi, DeltaEta, collision.centFT0C(), collision.posZ(), T_Sign, associate.sign(), 1. / TrackEff(hEffProtons, associate.sign(), associate.pt()));
                }
              }
            }
          }
        }
      }
    }
    // End of the V0-Track Correlations
  }

  void processMixed(MyFilteredCollisions const&, MyFilteredV0s const&, MyFilteredTracks const&)
  {

    // Start of the Mixed-events Correlations
    for (const auto& [coll_1, v0_1, coll_2, track_2] : pair) {
      for (const auto& [trigger, associate] : soa::combinations(soa::CombinationsFullIndexPolicy(v0_1, track_2))) {
        if (V0Filters(trigger) && TrackFilters(trigger, associate)) {
          if (FakeV0Filter(trigger, associate)) {

            T_Sign = V0Sign(trigger);
            if (T_Sign == 1) {
              CandMass = trigger.mLambda();
            } else if (T_Sign == -1) {
              CandMass = trigger.mAntiLambda();
            }

            A_PID = TrackPID(associate);
            DeltaPhi = DeltaPhiShift(trigger.phi(), associate.phi());
            DeltaEta = trigger.eta() - associate.eta();

            if (CandMass >= massLambda - 4 * DGaussSigma && CandMass <= massLambda + 4 * DGaussSigma) {
              if (A_PID[0] == 0.0) { // Pions
                MECorrRegistry.fill(HIST("hMixLambdaPion_SGNL"), DeltaPhi, DeltaEta, coll_1.centFT0C(), coll_1.posZ(), T_Sign, associate.sign(), 1. / TrackEff(hEffPions, associate.sign(), associate.pt()));
              } else if (A_PID[0] == 1.0) { // Kaons
                MECorrRegistry.fill(HIST("hMixLambdaKaon_SGNL"), DeltaPhi, DeltaEta, coll_1.centFT0C(), coll_1.posZ(), T_Sign, associate.sign(), 1. / TrackEff(hEffKaons, associate.sign(), associate.pt()));
              } else if (A_PID[0] == 2.0) { // Protons
                MECorrRegistry.fill(HIST("hMixLambdaProton_SGNL"), DeltaPhi, DeltaEta, coll_1.centFT0C(), coll_1.posZ(), T_Sign, associate.sign(), 1. / TrackEff(hEffProtons, associate.sign(), associate.pt()));
              }
            } else if (CandMass >= massLambda - 8 * DGaussSigma && CandMass <= massLambda + 8 * DGaussSigma) {
              if (A_PID[0] == 0.0) { // Pions
                MECorrRegistry.fill(HIST("hMixLambdaPion_SB"), DeltaPhi, DeltaEta, coll_1.centFT0C(), coll_1.posZ(), T_Sign, associate.sign(), 1. / TrackEff(hEffPions, associate.sign(), associate.pt()));
              } else if (A_PID[0] == 1.0) { // Kaons
                MECorrRegistry.fill(HIST("hMixLambdaKaon_SB"), DeltaPhi, DeltaEta, coll_1.centFT0C(), coll_1.posZ(), T_Sign, associate.sign(), 1. / TrackEff(hEffKaons, associate.sign(), associate.pt()));
              } else if (A_PID[0] == 2.0) { // Protons
                MECorrRegistry.fill(HIST("hMixLambdaProton_SB"), DeltaPhi, DeltaEta, coll_1.centFT0C(), coll_1.posZ(), T_Sign, associate.sign(), 1. / TrackEff(hEffProtons, associate.sign(), associate.pt()));
              }
            }
          }
        }
      }
    }
    // End of the Mixed-events Correlations
  }

  void processMCGen(MyFilteredMCGenCollision const&, MyFilteredMCParticles const& particles)
  {

    // Start of the Monte-Carlo generated QA
    for (const auto& particle : particles) {
      if (particle.isPhysicalPrimary()) {

        if (particle.pdgCode() == kPiPlus) { // Pos pions
          MCRegistry.fill(HIST("hGenPionP"), particle.pt());
        } else if (particle.pdgCode() == kPiMinus) { // Neg pions
          MCRegistry.fill(HIST("hGenPionN"), particle.pt());
        } else if (particle.pdgCode() == kKPlus) { // Pos kaons
          MCRegistry.fill(HIST("hGenKaonP"), particle.pt());
        } else if (particle.pdgCode() == kKMinus) { // Neg kaons
          MCRegistry.fill(HIST("hGenKaonN"), particle.pt());
        } else if (particle.pdgCode() == kProton) { // Pos protons
          MCRegistry.fill(HIST("hGenProtonP"), particle.pt());
        } else if (particle.pdgCode() == kProtonBar) { // Neg protons
          MCRegistry.fill(HIST("hGenProtonN"), particle.pt());
        }
      }
    }
    // End of the Monte-Carlo generated QA
  }

  void processMCRec(MyFilteredMCRecCollision const& collision, MyFilteredMCTracks const& tracks, aod::McCollisions const&, aod::McParticles const&)
  {

    if (!collision.has_mcCollision()) {
      return;
    }

    // Start of the Monte-Carlo reconstructed QA
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }

      auto particle = track.mcParticle();
      if (particle.isPhysicalPrimary()) {

        if (particle.pdgCode() == kPiPlus) { // Pos pions
          MCRegistry.fill(HIST("hRecPionP"), track.pt());
        } else if (particle.pdgCode() == kPiMinus) { // Neg pions
          MCRegistry.fill(HIST("hRecPionN"), track.pt());
        } else if (particle.pdgCode() == kKPlus) { // Pos kaons
          MCRegistry.fill(HIST("hRecKaonP"), track.pt());
        } else if (particle.pdgCode() == kKMinus) { // Neg kaons
          MCRegistry.fill(HIST("hRecKaonN"), track.pt());
        } else if (particle.pdgCode() == kProton) { // Pos protons
          MCRegistry.fill(HIST("hRecProtonP"), track.pt());
        } else if (particle.pdgCode() == kProtonBar) { // Neg protons
          MCRegistry.fill(HIST("hRecProtonN"), track.pt());
        }
      }
    }
    // End of the Monte-Carlo reconstructed QA
  }

  PROCESS_SWITCH(ThreePartCorr, processSame, "Process same-event correlations", true);
  PROCESS_SWITCH(ThreePartCorr, processMixed, "Process mixed-event correlations", true);
  PROCESS_SWITCH(ThreePartCorr, processMCGen, "Process Monte-Carlo, generator level", false);
  PROCESS_SWITCH(ThreePartCorr, processMCRec, "Process Monte-Carlo, reconstructed level", false);

  //================================================================================================================================================================================================================

  Double_t DeltaPhiShift(Double_t TriggerPhi, Double_t AssociatePhi)
  {

    Double_t dPhi = TriggerPhi - AssociatePhi;

    if (dPhi < (-1. / 2) * M_PI) {
      dPhi = dPhi + 2 * M_PI;
    } else if (dPhi > (3. / 2) * M_PI) {
      dPhi = dPhi - 2 * M_PI;
    }

    return dPhi;
  }

  Double_t TrackEff(TH1D** Efficiencies, Int_t Sign, Double_t pT)
  {

    Int_t Index = -999;
    if (Sign > 0) {
      Index = 0;
    } else if (Sign < 0) {
      Index = 1;
    }

    return Efficiencies[Index]->GetBinContent(Efficiencies[Index]->FindBin(pT));
  }

  template <class V0Cand>
  Int_t V0Sign(const V0Cand& V0)
  {

    if (TMath::Abs(V0.mLambda() - massLambda) <= TMath::Abs(V0.mAntiLambda() - massLambda)) {
      return 1;
    } else if (TMath::Abs(V0.mLambda() - massLambda) > TMath::Abs(V0.mAntiLambda() - massLambda)) {
      return -1;
    }

    return 0;
  }

  template <class TrackCand>
  Double_t* TrackPID(const TrackCand& Track)
  {

    static Double_t ID[2]; // {PID, NSigma}

    Double_t NSigmaTPC[3], NSigma[3];
    Double_t NSigmaTOF[3] = {0.0, 0.0, 0.0};
    NSigmaTPC[0] = Track.tpcNSigmaPi();
    NSigmaTPC[1] = Track.tpcNSigmaKa();
    NSigmaTPC[2] = Track.tpcNSigmaPr();
    if (Track.hasTOF()) {
      NSigmaTOF[0] = Track.tofNSigmaPi();
      NSigmaTOF[1] = Track.tofNSigmaKa();
      NSigmaTOF[2] = Track.tofNSigmaPr();
    }

    NSigma[0] = TMath::Sqrt(pow(NSigmaTPC[0], 2) + pow(NSigmaTOF[0], 2));
    NSigma[1] = TMath::Sqrt(pow(NSigmaTPC[1], 2) + pow(NSigmaTOF[1], 2));
    NSigma[2] = TMath::Sqrt(pow(NSigmaTPC[2], 2) + pow(NSigmaTOF[2], 2));

    if (NSigma[0] <= std::min(NSigma[1], NSigma[2])) { // Pions
      ID[0] = 0.0;
      ID[1] = NSigma[0];
    } else if (NSigma[1] <= std::min(NSigma[0], NSigma[2])) { // Kaons
      ID[0] = 1.0;
      ID[1] = NSigma[1];
    } else if (NSigma[2] <= std::min(NSigma[0], NSigma[1])) { // Protons
      ID[0] = 2.0;
      ID[1] = NSigma[2];
    }

    return ID;
  }

  template <class V0Cand>
  Bool_t V0Filters(const V0Cand& V0)
  {

    if (V0Sign(V0) == 1) {
      const auto& posDaughter = V0.template posTrack_as<MyFilteredTracks>();
      if (TMath::Abs(posDaughter.tpcNSigmaPr()) > 4.0) {
        return kFALSE;
      }
    } else if (V0Sign(V0) == -1) {
      const auto& negDaughter = V0.template negTrack_as<MyFilteredTracks>();
      if (TMath::Abs(negDaughter.tpcNSigmaPr()) > 4.0) {
        return kFALSE;
      }
    }

    return kTRUE;
  }

  template <class V0Cand, class TrackCand>
  Bool_t TrackFilters(const V0Cand& V0, const TrackCand& Track)
  {

    if (Track.globalIndex() == V0.posTrackId() || Track.globalIndex() == V0.negTrackId()) {
      return kFALSE;
    }
    if (TrackPID(Track)[1] > 4.0) {
      return kFALSE;
    }

    return kTRUE;
  }

  template <class V0Cand, class TrackCand>
  Bool_t FakeV0Filter(const V0Cand& V0, const TrackCand& Track)
  {

    if (ConfFilterSwitch) {

      if (TrackPID(Track)[0] == 1.0) { // Kaons
        return kTRUE;
      }

      std::array<float, 2> MassArray;
      std::array<float, 3> DMomArray;
      std::array<float, 3> AMomArray = Track.pVector();
      if (TrackPID(Track)[0] == 0.0) {
        MassArray = {o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged};

        if (V0Sign(V0) == 1 && Track.sign() == -1) { // Lambda - Pi_min
          const auto& dTrack = V0.template posTrack_as<MyFilteredTracks>();
          DMomArray = dTrack.pVector();
        } else if (V0Sign(V0) == -1 && Track.sign() == 1) { // Antilambda - Pi_plus
          const auto& dTrack = V0.template negTrack_as<MyFilteredTracks>();
          DMomArray = dTrack.pVector();
        }
      } else if (TrackPID(Track)[0] == 2.0) {
        MassArray = {o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton};

        if (V0Sign(V0) == 1 && Track.sign() == 1) { // Lambda - Proton
          const auto& dTrack = V0.template negTrack_as<MyFilteredTracks>();
          DMomArray = dTrack.pVector();
        } else if (V0Sign(V0) == -1 && Track.sign() == -1) { // Antilambda - Antiproton
          const auto& dTrack = V0.template posTrack_as<MyFilteredTracks>();
          DMomArray = dTrack.pVector();
        }
      }

      Double_t M = RecoDecay::m(std::array{DMomArray, AMomArray}, MassArray);
      if (M >= massLambda - 4 * DGaussSigma && M <= massLambda + 4 * DGaussSigma) {
        return kFALSE;
      }
    }

    return kTRUE;
  }
};

//==================================================================================================================================================================================================================

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<ThreePartCorr>(cfgc)};
  return workflow;
}

//==================================================================================================================================================================================================================
