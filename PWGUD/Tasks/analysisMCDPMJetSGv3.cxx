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
/// \file AnalysisMCDPMJetSGv3.cxx
/// \brief Process MC DPMJet events for inclusive studies.
///
/// \author Simone Ragoni

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
// #include "TDatabasePDG.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/DataModel/UDTables.h"
// #include "TLorentzVector.h"
// #include "TVector3.h"
#include "Math/LorentzVector.h" // ROOT::Math::LorentzVector
#include "Math/PxPyPzM4D.h"     // ROOT::Math::PxPyPzM4D
#include "TMath.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct AnalysisMCDPMJetSGv3 {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // TDatabasePDG* fPDG = TDatabasePDG::Instance();

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histo"};

  // using myCompleteTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels>;
  // using myFilteredTracks = soa::Filtered<myCompleteTracks>;
  using BCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  Preslice<o2::aod::McParticles> perMcCollision = o2::aod::mcparticle::mcCollisionId;
  // using MCTCs = soa::Join<aod::Tracks, aod::TracksExtra, /*aod::TracksCov,*/ aod::TracksDCA, aod::TrackSelection,
  //                         aod::McTrackLabels,
  //                         aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
  //                         aod::TOFSignal, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  // using MCTC = MCTCs::iterator;
  // define abbreviations
  using CCs = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDMcCollsLabels>;
  using CC = CCs::iterator;
  using MCparticles = aod::UDMcParticles::iterator;
  using TCs = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA, aod::UDMcTrackLabels>;
  // using TCs = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksPID, aod::UDMcTrackLabels>;
  using TC = TCs::iterator;
  using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>;

  double massPion = 0.;
  double massKaon = 0.;
  double massProton = 0.;
  const int codePion = 211;
  const int codeKaon = 321;
  const int codeProton = 2212;

  void init(InitContext const&)
  {
    // TParticlePDG* pionPDG = fPDG->GetParticle(codePion);
    // if (pionPDG != nullptr) {
    //   massPion = pionPDG->Mass();
    // }
    // TParticlePDG* kaonPDG = fPDG->GetParticle(codeKaon);
    // if (kaonPDG != nullptr) {
    //   massKaon = kaonPDG->Mass();
    // }
    // TParticlePDG* protonPDG = fPDG->GetParticle(codeProton);
    // if (protonPDG != nullptr) {
    //   massProton = protonPDG->Mass();
    // }
    massPion = o2::constants::physics::MassPionCharged;
    massKaon = o2::constants::physics::MassKaonCharged;
    massProton = o2::constants::physics::MassProton;

    // define axes you want to use
    const AxisSpec axisCounter{10, 0, 10, ""};
    const AxisSpec axisEta{100, -1.5, +1.5, "#eta"};
    const AxisSpec axisPt{5000, 0, 5, "p_{T}"};
    const AxisSpec axisPtSmall{1000, 0, 1, "p_{T}"};
    const AxisSpec axisMass{nBinsPt, 0, 5, "m_{#pi#pi}"};
    const AxisSpec axisMassSmall{nBinsPt, 0, 2, "m_{#pi#pi}"};
    const AxisSpec axisDeltaPt{100, -1.0, +1.0, "#Delta(p_{T})"};
    const AxisSpec axisBC{1000, -10000.0, +10000.0, "BCs"};
    const AxisSpec axisBCext{100000, -10000000.0, +10000000.0, "BCs"};
    const AxisSpec axisCosTheta{100, -1.0, +1.0, "cos#theta"};
    const AxisSpec axisPhi{600, -o2::constants::math::PI, -o2::constants::math::PI, "#varphi"};

    // create histograms
    histos.add("hdEdx", "p vs dE/dx Signal", kTH2F, {{100, 0.0, 3.0}, {1000, 0.0, 2000.0}});
    histos.add("hSigmaPion", "p vs dE/dx sigma pion   TPC      ", kTH2F, {{300, 0.0, 3.0}, {200, -10.0, 10.0}});
    histos.add("hSigmaPionTruth", "p vs dE/dx sigma pion   TPC truth", kTH2F, {{300, 0.0, 3.0}, {200, -10.0, 10.0}});
    histos.add("hSigmaPionTOF", "p vs dE/dx sigma pion   TOF      ", kTH2F, {{300, 0.0, 3.0}, {200, -10.0, 10.0}});
    histos.add("hSigmaPionTruthTOF", "p vs dE/dx sigma pion   TOF truth", kTH2F, {{300, 0.0, 3.0}, {200, -10.0, 10.0}});
    histos.add("hSigmaKaon", "p vs dE/dx sigma kaon   TPC      ", kTH2F, {{300, 0.0, 3.0}, {200, -10.0, 10.0}});
    histos.add("hSigmaKaonTruth", "p vs dE/dx sigma kaon   TPC truth", kTH2F, {{300, 0.0, 3.0}, {200, -10.0, 10.0}});
    histos.add("hSigmaKaonTOF", "p vs dE/dx sigma kaon   TOF      ", kTH2F, {{300, 0.0, 3.0}, {200, -10.0, 10.0}});
    histos.add("hSigmaKaonTruthTOF", "p vs dE/dx sigma kaon   TOF truth", kTH2F, {{300, 0.0, 3.0}, {200, -10.0, 10.0}});
    histos.add("hSigmaProton", "p vs dE/dx sigma proton TPC      ", kTH2F, {{300, 0.0, 3.0}, {200, -10.0, 10.0}});
    histos.add("hSigmaProtonTruth", "p vs dE/dx sigma proton TPC truth", kTH2F, {{300, 0.0, 3.0}, {200, -10.0, 10.0}});
    histos.add("hSigmaProtonTOF", "p vs dE/dx sigma proton TOF      ", kTH2F, {{300, 0.0, 3.0}, {200, -10.0, 10.0}});
    histos.add("hSigmaProtonTruthTOF", "p vs dE/dx sigma proton TOF truth", kTH2F, {{300, 0.0, 3.0}, {200, -10.0, 10.0}});

    histos.add("hVisibleMultiVsGeneratedMulti", "Multiplicity correlation", kTH2F, {{10000, -0.5, 9999.5}, {1000, -0.5, 999.5}});

    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    histos.add("ptResolution", "ptResolution", kTH2F, {axisPt, axisDeltaPt});

    histos.add("ptGeneratedPion", "ptGeneratedPion", kTH1F, {axisPt});
    histos.add("ptGeneratedKaon", "ptGeneratedKaon", kTH1F, {axisPt});
    histos.add("ptGeneratedProton", "ptGeneratedProton", kTH1F, {axisPt});
    histos.add("ptGeneratedPionAxE", "ptGeneratedPionAxE", kTH1F, {axisPt});
    histos.add("ptGeneratedKaonAxE", "ptGeneratedKaonAxE", kTH1F, {axisPt});
    histos.add("ptGeneratedProtonAxE", "ptGeneratedProtonAxE", kTH1F, {axisPt});
    histos.add("ptGeneratedProtonAxEPos", "ptGeneratedProtonAxEPos", kTH1F, {axisPt});
    histos.add("ptGeneratedProtonAxENeg", "ptGeneratedProtonAxENeg", kTH1F, {axisPt});
    histos.add("ptReconstructedPion", "ptReconstructedPion", kTH1F, {axisPt});
    histos.add("ptReconstructedKaon", "ptReconstructedKaon", kTH1F, {axisPt});
    histos.add("ptReconstructedProton", "ptReconstructedProton", kTH1F, {axisPt});
    histos.add("ptReconstructedProtonPos", "ptReconstructedProtonPos", kTH1F, {axisPt});
    histos.add("ptReconstructedProtonNeg", "ptReconstructedProtonNeg", kTH1F, {axisPt});
    histos.add("ptReconstructedPionTOF", "ptReconstructedPionTOF", kTH1F, {axisPt});
    histos.add("ptReconstructedKaonTOF", "ptReconstructedKaonTOF", kTH1F, {axisPt});
    histos.add("ptReconstructedProtonTOF", "ptReconstructedProtonTOF", kTH1F, {axisPt});

    histos.add("allreconstructedPFPion", "allreconstructedPFPion", kTH1F, {axisPt});
    histos.add("allreconstructedPFKaon", "allreconstructedPFKaon", kTH1F, {axisPt});
    histos.add("allreconstructedPFProton", "allreconstructedPFProton", kTH1F, {axisPt});
    histos.add("allreconstructedPFProtonPos", "allreconstructedPFProtonPos", kTH1F, {axisPt});
    histos.add("allreconstructedPFProtonNeg", "allreconstructedPFProtonNeg", kTH1F, {axisPt});
    histos.add("allreconstructedPFPionTOF", "allreconstructedPFPionTOF", kTH1F, {axisPt});
    histos.add("allreconstructedPFKaonTOF", "allreconstructedPFKaonTOF", kTH1F, {axisPt});
    histos.add("allreconstructedPFProtonTOF", "allreconstructedPFProtonTOF", kTH1F, {axisPt});

    histos.add("numberOfRecoCollisions", "numberOfRecoCollisions", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("numberOfRecoCollisions2", "numberOfRecoCollisions2", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("numberOfTracksMC", "numberOfTracksMC", kTH1F, {{100, -0.5f, 99.5f}});
    histos.add("numberOfTracksReco", "numberOfTracksReco", kTH1F, {{100, -0.5f, 99.5f}});

    histos.add("bcResolution", "bcResolution", kTH1F, {axisBC});
    histos.add("mcbcHistogram", "mcbcHistogram", kTH1F, {axisBCext});
    histos.add("bcHistogram", "bcHistogram", kTH1F, {axisBCext});
    histos.add("mcbcModuloOrbitHistogram", "mcbcModuloOrbitHistogram", kTH1F, {axisBC});
    histos.add("bcModuloOrbitHistogram", "bcModuloOrbitHistogram", kTH1F, {axisBC});
  }
  //-----------------------------------------------------------------------------------------------------------------------
  void processSim(aod::UDMcCollision const& mcCollision, aod::UDMcParticles const& mcParticles)
  {
    // histos.fill(HIST("eventCounter"), 0.5);

    // auto massPion = 0.;
    // TParticlePDG pionPDG = fPDG->GetParticle(codePion);
    // massPion = pionPDG.Mass();
    // auto massKaon = 0.;
    // TParticlePDG kaonPDG = fPDG->GetParticle(codeKaon);
    // massKaon = kaonPDG.Mass();
    // auto massProton = 0.;
    // TParticlePDG protonPDG = fPDG->GetParticle(codeProton);
    // massProton = protonPDG.Mass();
    histos.fill(HIST("numberOfTracksMC"), mcParticles.size());
    histos.fill(HIST("eventCounter"), mcCollision.size());
    // LOGF(info, "New event! mcParticles.size() = %d", mcParticles.size());

    int counterMC = 0;
    int counter = 0;
    for (const auto& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary())
        continue;
      counterMC += 1;
      // if(mcParticle.isPhysicalPrimary()) counterMC += 1;
      LorentzVectorM protoMC(
        mcParticle.px(),
        mcParticle.py(),
        mcParticle.pz(),
        massPion);
      double etaMax = 0.8;
      double ptMin = 0.1;
      if (std::fabs(protoMC.Eta()) < etaMax && protoMC.Pt() > ptMin) {
        counter += 1;
      }
      if (!mcParticle.isPhysicalPrimary())
        continue;
      // if(mcParticle.isPhysicalPrimary() && fabs(mcParticle.eta())<0.9){ // do this in the context of the MC loop ! (context matters!!!)
      // LorentzVectorM pMC;
      LorentzVectorM pMC(mcParticle.px(), mcParticle.py(), mcParticle.pz(), massPion);
      if (std::abs(mcParticle.pdgCode()) == codePion) {
        // histos.fill(HIST("ptGeneratedPion"), mcParticle.pt());
        // LorentzVectorM pMC(mcParticle.px(), mcParticle.py(), mcParticle.pz(), massPion);
        histos.fill(HIST("ptGeneratedPion"), pMC.Pt());
      }
      if (std::abs(mcParticle.pdgCode()) == codeKaon) {
        // histos.fill(HIST("ptGenerateKaon"), mcParticle.pt());
        // LorentzVectorM pMC(mcParticle.px(), mcParticle.py(), mcParticle.pz(), massKaon);
        pMC.SetM(massKaon);
        histos.fill(HIST("ptGeneratedKaon"), pMC.Pt());
      }
      if (std::abs(mcParticle.pdgCode()) == codeProton) {
        // histos.fill(HIST("ptGeneratedProton"), mcParticle.pt());
        // LorentzVectorM pMC(mcParticle.px(), mcParticle.py(), mcParticle.pz(), massProton);
        pMC.SetM(massProton);
        histos.fill(HIST("ptGeneratedProton"), pMC.Pt());
      }
      double yMax = 0.8;
      if (std::abs(pMC.Rapidity()) < yMax) {
        if (std::abs(mcParticle.pdgCode()) == codePion)
          histos.fill(HIST("ptGeneratedPionAxE"), pMC.Pt());
        if (std::abs(mcParticle.pdgCode()) == codeKaon)
          histos.fill(HIST("ptGeneratedKaonAxE"), pMC.Pt());
        if (std::abs(mcParticle.pdgCode()) == codeProton)
          histos.fill(HIST("ptGeneratedProtonAxE"), pMC.Pt());
        if (mcParticle.pdgCode() == codeProton) {
          histos.fill(HIST("ptGeneratedProtonAxEPos"), pMC.Pt());
        } else {
          histos.fill(HIST("ptGeneratedProtonAxENeg"), pMC.Pt());
        }
      }
    }
    histos.fill(HIST("hVisibleMultiVsGeneratedMulti"), counterMC, counter);
  }
  PROCESS_SWITCH(AnalysisMCDPMJetSGv3, processSim, "processSim", true);

  void processReco(CC const& collision,
                   TCs const& tracks,
                   // aod::UDMcCollisions const& /*mccollisions*/,
                   aod::UDMcParticles const& mcParticles)
  {
    histos.fill(HIST("numberOfRecoCollisions"), 88.);              // number of times coll was reco-ed
    histos.fill(HIST("numberOfRecoCollisions"), collision.size()); // number of times coll was reco-ed
    histos.fill(HIST("numberOfRecoCollisions2"), mcParticles.size());
    Partition<TCs> pvContributors = aod::udtrack::isPVContributor == true;
    pvContributors.bindTable(tracks);

    // auto massPion = 0.;
    // TParticlePDG pionPDG = fPDG->GetParticle(codePion);
    // massPion = pionPDG.Mass();
    // auto massKaon = 0.;
    // TParticlePDG kaonPDG = fPDG->GetParticle(codeKaon);
    // massKaon = kaonPDG.Mass();
    // auto massProton = 0.;
    // TParticlePDG protonPDG = fPDG->GetParticle(codeProton);
    // massProton = protonPDG.Mass();

    histos.fill(HIST("numberOfTracksReco"), tracks.size());
    // double etaMax = 0.8;
    double yMax = 0.8;
    double sigmaMax = 3.;
    double ptMin = 0.1;
    int nFindableMin = 70;
    double dcaZlimit = 2.;

    // int counter = 0;
    for (const auto& track : tracks) {
      if (track.isPVContributor()) {
        int nFindable = track.tpcNClsFindable();
        if (nFindable < nFindableMin) {
          continue;
        }
        // int NMinusFound = track.tpcNClsFindableMinusFound();
        // int NCluster = NFindable - NMinusFound;
        // if (NCluster < 70) {
        //   continue;
        // }
        if (track.pt() < ptMin) {
          continue;
        }
        if (!(std::abs(track.dcaZ()) < dcaZlimit)) {
          continue;
        }
        double dcaLimit = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
        if (!(std::abs(track.dcaXY()) < dcaLimit)) {
          continue;
        }

        double momentum = std::sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());
        double dEdx = track.tpcSignal();
        histos.fill(HIST("hdEdx"), momentum, dEdx);

        LorentzVectorM pion(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged);
        LorentzVectorM kaon(track.px(), track.py(), track.pz(), o2::constants::physics::MassKaonCharged);
        LorentzVectorM proton(track.px(), track.py(), track.pz(), o2::constants::physics::MassProton);
        auto nSigmaPi = -999.;
        auto nSigmaKa = -999.;
        auto nSigmaPr = -999.;
        auto nSigmaPiTOF = -999.;
        auto nSigmaKaTOF = -999.;
        auto nSigmaPrTOF = -999.;
        // This section makes templates
        if (track.hasTPC()) {
          nSigmaPi = track.tpcNSigmaPi();
          nSigmaKa = track.tpcNSigmaKa();
          nSigmaPr = track.tpcNSigmaPr();
          if (std::abs(nSigmaPi) < sigmaMax && std::abs(pion.Rapidity()) < yMax) {
            histos.fill(HIST("hSigmaPion"), track.pt(), nSigmaPi);
            if (track.has_udMcParticle()) {
              auto mcParticle = track.udMcParticle();
              // if(abs(mcParticle.pdgCode())==codePion && mcParticle.isPhysicalPrimary()) howManyPionsHavePionMCandPrimaries += 1;
              if (std::abs(mcParticle.pdgCode()) == codePion) {
                histos.fill(HIST("hSigmaPionTruth"), track.pt(), nSigmaPi);
                histos.fill(HIST("allreconstructedPFPion"), track.pt());
                if (mcParticle.isPhysicalPrimary()) {
                  histos.fill(HIST("ptReconstructedPion"), track.pt());
                }
              }
            }
          }
          if (std::abs(nSigmaKa) < sigmaMax && std::abs(kaon.Rapidity()) < yMax) {
            histos.fill(HIST("hSigmaKaon"), track.pt(), nSigmaKa);
            if (track.has_udMcParticle()) {
              auto mcParticle = track.udMcParticle();
              if (std::abs(mcParticle.pdgCode()) == codeKaon) {
                histos.fill(HIST("hSigmaKaonTruth"), track.pt(), nSigmaKa);
                histos.fill(HIST("allreconstructedPFKaon"), track.pt());
                if (mcParticle.isPhysicalPrimary()) {
                  histos.fill(HIST("ptReconstructedKaon"), track.pt());
                }
              }
            }
          }
          if (std::abs(nSigmaPr) < sigmaMax && std::abs(proton.Rapidity()) < yMax) {
            histos.fill(HIST("hSigmaProton"), track.pt(), nSigmaPr);
            if (track.has_udMcParticle()) {
              auto mcParticle = track.udMcParticle();
              if (std::abs(mcParticle.pdgCode()) == codeProton) {
                histos.fill(HIST("hSigmaProtonTruth"), track.pt(), nSigmaPr);
                histos.fill(HIST("allreconstructedPFProton"), track.pt());
                if (mcParticle.pdgCode() == codeProton) {
                  histos.fill(HIST("allreconstructedPFProtonPos"), track.pt());
                } else {
                  histos.fill(HIST("allreconstructedPFProtonNeg"), track.pt());
                }
                if (mcParticle.isPhysicalPrimary()) {
                  histos.fill(HIST("ptReconstructedProton"), track.pt());
                  if (mcParticle.pdgCode() == codeProton) {
                    histos.fill(HIST("ptReconstructedProtonPos"), track.pt());
                  } else {
                    histos.fill(HIST("ptReconstructedProtonNeg"), track.pt());
                  }
                }
              }
            }
          }
        }
        if (track.hasTPC() && track.hasTOF()) {
          // if (track.hasTOF()) {
          nSigmaPiTOF = track.tofNSigmaPi();
          nSigmaKaTOF = track.tofNSigmaKa();
          nSigmaPrTOF = track.tofNSigmaPr();
          if (std::abs(nSigmaPiTOF) < sigmaMax && std::abs(pion.Rapidity()) < yMax) {
            histos.fill(HIST("hSigmaPionTOF"), track.pt(), nSigmaPiTOF);
            if (track.has_udMcParticle()) {
              auto mcParticle = track.udMcParticle();
              // if(abs(mcParticle.pdgCode())==codePion && mcParticle.isPhysicalPrimary()) howManyPionsHavePionMCandPrimaries += 1;
              if (std::abs(mcParticle.pdgCode()) == codePion) {
                histos.fill(HIST("hSigmaPionTruthTOF"), track.pt(), nSigmaPiTOF);
                histos.fill(HIST("allreconstructedPFPionTOF"), track.pt());
                if (mcParticle.isPhysicalPrimary()) {
                  histos.fill(HIST("ptReconstructedPionTOF"), track.pt());
                }
              }
            }
          }
          if (std::abs(nSigmaKaTOF) < sigmaMax && std::abs(kaon.Rapidity()) < yMax) {
            histos.fill(HIST("hSigmaKaonTOF"), track.pt(), nSigmaKaTOF);
            if (track.has_udMcParticle()) {
              auto mcParticle = track.udMcParticle();
              if (std::abs(mcParticle.pdgCode()) == codeKaon) {
                histos.fill(HIST("hSigmaKaonTruthTOF"), track.pt(), nSigmaKaTOF);
                histos.fill(HIST("allreconstructedPFKaonTOF"), track.pt());
                if (mcParticle.isPhysicalPrimary()) {
                  histos.fill(HIST("ptReconstructedKaonTOF"), track.pt());
                }
              }
            }
          }
          if (std::abs(nSigmaPrTOF) < sigmaMax && std::abs(proton.Rapidity()) < yMax) {
            histos.fill(HIST("hSigmaProtonTOF"), track.pt(), nSigmaPrTOF);
            if (track.has_udMcParticle()) {
              auto mcParticle = track.udMcParticle();
              if (std::abs(mcParticle.pdgCode()) == codeProton) {
                histos.fill(HIST("hSigmaProtonTruthTOF"), track.pt(), nSigmaPrTOF);
                histos.fill(HIST("allreconstructedPFProtonTOF"), track.pt());
                if (mcParticle.isPhysicalPrimary()) {
                  histos.fill(HIST("ptReconstructedProtonTOF"), track.pt());
                }
              }
            }
          }
        }
        // counter++;
        // histos.fill(HIST("hVisibleMultiVsGeneratedMulti"), counterMC, counter);
        // histos.fill(HIST("hVisibleMultiVsGeneratedMulti"), mcParticles.size(), counter);
      }
    } // track loop

    // } // collision loop
  }
  PROCESS_SWITCH(AnalysisMCDPMJetSGv3, processReco, "processReco", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AnalysisMCDPMJetSGv3>(cfgc)};
}
