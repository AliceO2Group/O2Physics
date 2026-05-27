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
/// \file ptSpectraInclusiveUpc.cxx
/// \executable o2-analysis-ud-pt-spectra-inclusive-upc
/// \brief Task for the of pT spectra of pions, kaons and protons in inclusive UPC events.
///        Used to obtain the templates for the DCA_xy fits for the primary fractions
///
/// \author Andrea Giovanni Riffero andrea.giovanni.riffero@cern.ch

#include "PWGUD/DataModel/UDTables.h"

#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PxPyPzM4D.h>
#include <TMCProcess.h>
#include <TPDGCode.h>

#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PtSpectraInclusiveUpc {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pT histos"};
  Configurable<int> nBinsDCAxy{"nBinsDCAxy", 100, "N bins in DCA_{xy} histos"};
  Configurable<bool> applyKineCutsInGen{"applyKineCutsInGen", false, "Apply kinematic cuts in the generated level"};

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  Preslice<o2::aod::McParticles> perMcCollision = o2::aod::mcparticle::mcCollisionId;
  // define abbreviations
  using CCs = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDMcCollsLabels>;
  using CC = CCs::iterator;
  using TCs = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA, aod::UDMcTrackLabels>;
  using TC = TCs::iterator;
  using LorentzVectorM = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>;

  const double etaMax = 0.9;
  const double yMax = 0.9;
  const double ptMin = 0.1;
  const int nFindableMin = 70;
  const double sigmaMax = 3.;
  const double dcaZlimit = 2.;

  void init(InitContext const&)
  {
    // axes
    const AxisSpec axisPt{nBinsPt, 0, 5, "#it{p}_{T} GeV/#it{c}"};
    const AxisSpec axisEventCounter{2, 0.5, 2.5, "Event type"};
    const AxisSpec axisDCAxy{nBinsDCAxy, -0.6, 0.6, "DCA_{xy} cm"};

    // histograms
    histos.add("ptGeneratedPion", "ptGeneratedPion", kTH1F, {axisPt});
    histos.add("ptGeneratedKaon", "ptGeneratedKaon", kTH1F, {axisPt});
    histos.add("ptGeneratedProton", "ptGeneratedProton", kTH1F, {axisPt});

    histos.add("ptReconstructedTPCPion", "ptReconstructedTPCPion", kTH1F, {axisPt});
    histos.add("ptReconstructedTPCKaon", "ptReconstructedTPCKaon", kTH1F, {axisPt});
    histos.add("ptReconstructedTPCProton", "ptReconstructedTPCProton", kTH1F, {axisPt});

    histos.add("ptReconstructedTOFPion", "ptReconstructedTOFPion", kTH1F, {axisPt});
    histos.add("ptReconstructedTOFKaon", "ptReconstructedTOFKaon", kTH1F, {axisPt});
    histos.add("ptReconstructedTOFProton", "ptReconstructedTOFProton", kTH1F, {axisPt});

    histos.add("ptDataTPCPion", "ptDataTPCPion", kTH1F, {axisPt});
    histos.add("ptDataTPCKaon", "ptDataTPCKaon", kTH1F, {axisPt});
    histos.add("ptDataTPCProton", "ptDataTPCProton", kTH1F, {axisPt});

    histos.add("ptDataTOFPion", "ptDataTOFPion", kTH1F, {axisPt});
    histos.add("ptDataTOFKaon", "ptDataTOFKaon", kTH1F, {axisPt});
    histos.add("ptDataTOFProton", "ptDataTOFProton", kTH1F, {axisPt});

    histos.add("myEventCounter", "myEventCounter", kTH1F, {axisEventCounter});
    histos.add("DCAxy_primary_pions", "DCAxy_primary_pions", kTH1F, {axisDCAxy});
    histos.add("DCAxy_secondary_pions", "DCAxy_secondary_pions", kTH1F, {axisDCAxy});
    histos.add("DCAxy_primary_kaons", "DCAxy_primary_kaons", kTH1F, {axisDCAxy});
    histos.add("DCAxy_secondary_kaons", "DCAxy_secondary_kaons", kTH1F, {axisDCAxy});
    histos.add("DCAxy_primary_protons", "DCAxy_primary_protons", kTH1F, {axisDCAxy});
    histos.add("DCAxy_secondary_protons", "DCAxy_secondary_protons", kTH1F, {axisDCAxy});
    histos.add("DCAxy_material_protons", "DCAxy_material_protons", kTH1F, {axisDCAxy});
    histos.add("DCAxy_data_pions", "DCAxy_data_pions", kTH1F, {axisDCAxy});
    histos.add("DCAxy_data_kaons", "DCAxy_data_kaons", kTH1F, {axisDCAxy});
    histos.add("DCAxy_data_protons", "DCAxy_data_protons", kTH1F, {axisDCAxy});
  }

  void processSim(aod::UDMcCollision const&, aod::UDMcParticles const& mcParticles)
  {

    for (const auto& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary())
        continue;

      LorentzVectorM pMC(mcParticle.px(), mcParticle.py(), mcParticle.pz(), o2::constants::physics::MassPionCharged);

      if (applyKineCutsInGen) {
        if (std::fabs(pMC.Eta()) > etaMax)
          continue;

        if (std::fabs(pMC.Rapidity()) > yMax)
          continue;

        if (pMC.Pt() < ptMin)
          continue;
      }

      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus) {
        histos.fill(HIST("ptGeneratedPion"), pMC.Pt());
      }

      if (std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus) {
        pMC.SetM(o2::constants::physics::MassKaonCharged);
        histos.fill(HIST("ptGeneratedKaon"), pMC.Pt());
      }

      if (std::abs(mcParticle.pdgCode()) == PDG_t::kProton) {
        pMC.SetM(o2::constants::physics::MassProton);
        histos.fill(HIST("ptGeneratedProton"), pMC.Pt());
      }

      histos.fill(HIST("myEventCounter"), 1); // gen event
    }
  }

  void processReco(CC const&, TCs const& tracks, aod::UDMcParticles const&)
  {

    double dcaXyLimit = 0;

    auto nSigmaPi = -999.;
    auto nSigmaKa = -999.;
    auto nSigmaPr = -999.;

    LorentzVectorM* pion = new LorentzVectorM();
    LorentzVectorM* kaon = new LorentzVectorM();
    LorentzVectorM* proton = new LorentzVectorM();

    for (const auto& track : tracks) {
      if (!track.isPVContributor()) {
        continue;
      }

      if (track.tpcNClsFindable() < nFindableMin) {
        continue;
      }

      if (track.pt() < ptMin) {
        continue;
      }

      if (!(std::abs(track.dcaZ()) < dcaZlimit)) {
        continue;
      }

      dcaXyLimit = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
      if (!(std::abs(track.dcaXY()) < dcaXyLimit)) {
        continue;
      }

      pion->SetPx(track.px());
      pion->SetPy(track.py());
      pion->SetPz(track.pz());
      pion->SetM(o2::constants::physics::MassPionCharged);

      kaon->SetPx(track.px());
      kaon->SetPy(track.py());
      kaon->SetPz(track.pz());
      kaon->SetM(o2::constants::physics::MassKaonCharged);

      proton->SetPx(track.px());
      proton->SetPy(track.py());
      proton->SetPz(track.pz());
      proton->SetM(o2::constants::physics::MassProton);

      if (!track.has_udMcParticle()) {
        continue;
      }
      auto mcParticle = track.udMcParticle();

      bool hasTpc = false;
      // TPC tracks
      if (track.hasTPC()) {
        hasTpc = true;
        nSigmaPi = track.tpcNSigmaPi();
        nSigmaKa = track.tpcNSigmaKa();
        nSigmaPr = track.tpcNSigmaPr();

        if (std::abs(nSigmaPi) < sigmaMax) {
          if (std::abs(pion->Rapidity()) > yMax) {
            continue;
          }

          if (mcParticle.isPhysicalPrimary()) {
            histos.fill(HIST("ptReconstructedTPCPion"), pion->Pt());
            histos.fill(HIST("DCAxy_primary_pions"), track.dcaXY());
          } else {
            histos.fill(HIST("DCAxy_secondary_pions"), track.dcaXY());
          }
        }
        if (std::abs(nSigmaKa) < sigmaMax) {
          if (std::abs(kaon->Rapidity()) > yMax) {
            continue;
          }

          if (mcParticle.isPhysicalPrimary()) {
            histos.fill(HIST("ptReconstructedTPCKaon"), kaon->Pt());
            histos.fill(HIST("DCAxy_primary_kaons"), track.dcaXY());
          } else {
            histos.fill(HIST("DCAxy_secondary_kaons"), track.dcaXY());
          }
        }

        if (std::abs(nSigmaPr) < sigmaMax) {
          if (std::abs(proton->Rapidity()) > yMax) {
            continue;
          }

          if (mcParticle.isPhysicalPrimary()) {
            histos.fill(HIST("ptReconstructedTPCProton"), proton->Pt());
            histos.fill(HIST("DCAxy_primary_protons"), track.dcaXY());
          } else {
            if (mcParticle.getProcess() == kPDecay) {
              histos.fill(HIST("DCAxy_secondary_protons"), track.dcaXY());
            } else {
              histos.fill(HIST("DCAxy_material_protons"), track.dcaXY());
            }
          }
        }
      }

      // TPC tracks
      if (track.hasTOF()) {
        nSigmaPi = track.tofNSigmaPi();
        nSigmaKa = track.tofNSigmaKa();
        nSigmaPr = track.tofNSigmaPr();

        if (std::abs(nSigmaPi) < sigmaMax) {
          if (std::abs(pion->Rapidity()) > yMax) {
            continue;
          }

          if (mcParticle.isPhysicalPrimary()) {
            histos.fill(HIST("ptReconstructedTOFPion"), pion->Pt());
            if (!hasTpc)
              histos.fill(HIST("DCAxy_primary_pions"), track.dcaXY());
          } else {
            if (!hasTpc)
              histos.fill(HIST("DCAxy_secondary_pions"), track.dcaXY());
          }
        }
        if (std::abs(nSigmaKa) < sigmaMax) {
          if (std::abs(kaon->Rapidity()) > yMax) {
            continue;
          }

          if (mcParticle.isPhysicalPrimary()) {
            histos.fill(HIST("ptReconstructedTOFKaon"), kaon->Pt());
            if (!hasTpc)
              histos.fill(HIST("DCAxy_primary_kaons"), track.dcaXY());
          } else {
            if (!hasTpc)
              histos.fill(HIST("DCAxy_secondary_kaons"), track.dcaXY());
          }
        }

        if (std::abs(nSigmaPr) < sigmaMax) {
          if (std::abs(proton->Rapidity()) > yMax) {
            continue;
          }

          if (mcParticle.isPhysicalPrimary()) {
            histos.fill(HIST("ptReconstructedTOFProton"), proton->Pt());
            if (!hasTpc)
              histos.fill(HIST("DCAxy_primary_protons"), track.dcaXY());
          } else {
            if (!hasTpc) {
              if (mcParticle.getProcess() == kPDecay) {
                histos.fill(HIST("DCAxy_secondary_protons"), track.dcaXY());
              } else {
                histos.fill(HIST("DCAxy_material_protons"), track.dcaXY());
              }
            }
          }
        }
      }
    }
  }

  void processData(CC const&, TCs const& tracks)
  {

    double dcaXyLimit = 0;

    auto nSigmaPi = -999.;
    auto nSigmaKa = -999.;
    auto nSigmaPr = -999.;

    LorentzVectorM* pion = new LorentzVectorM();
    LorentzVectorM* kaon = new LorentzVectorM();
    LorentzVectorM* proton = new LorentzVectorM();

    for (const auto& track : tracks) {
      if (!track.isPVContributor()) {
        continue;
      }

      if (track.tpcNClsFindable() < nFindableMin) {
        continue;
      }

      if (track.pt() < ptMin) {
        continue;
      }

      if (!(std::abs(track.dcaZ()) < dcaZlimit)) {
        continue;
      }

      dcaXyLimit = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
      if (!(std::abs(track.dcaXY()) < dcaXyLimit)) {
        continue;
      }

      pion->SetPx(track.px());
      pion->SetPy(track.py());
      pion->SetPz(track.pz());
      pion->SetM(o2::constants::physics::MassPionCharged);

      kaon->SetPx(track.px());
      kaon->SetPy(track.py());
      kaon->SetPz(track.pz());
      kaon->SetM(o2::constants::physics::MassKaonCharged);

      proton->SetPx(track.px());
      proton->SetPy(track.py());
      proton->SetPz(track.pz());
      proton->SetM(o2::constants::physics::MassProton);

      bool hasTpc = false;
      // TPC tracks
      if (track.hasTPC()) {
        hasTpc = true;
        nSigmaPi = track.tpcNSigmaPi();
        nSigmaKa = track.tpcNSigmaKa();
        nSigmaPr = track.tpcNSigmaPr();

        if (std::abs(nSigmaPi) < sigmaMax) {
          if (std::abs(pion->Rapidity()) > yMax) {
            continue;
          }
          histos.fill(HIST("ptDataTPCPion"), pion->Pt());
          histos.fill(HIST("DCAxy_data_pions"), track.dcaXY());
        }

        if (std::abs(nSigmaKa) < sigmaMax) {
          if (std::abs(kaon->Rapidity()) > yMax) {
            continue;
          }
          histos.fill(HIST("ptDataTPCKaon"), kaon->Pt());
          histos.fill(HIST("DCAxy_data_kaons"), track.dcaXY());
        }

        if (std::abs(nSigmaPr) < sigmaMax) {
          if (std::abs(proton->Rapidity()) > yMax) {
            continue;
          }
          histos.fill(HIST("ptDataTPCProton"), proton->Pt());
          histos.fill(HIST("DCAxy_data_protons"), track.dcaXY());
        }
      }

      // TPC tracks
      if (track.hasTOF()) {
        nSigmaPi = track.tofNSigmaPi();
        nSigmaKa = track.tofNSigmaKa();
        nSigmaPr = track.tofNSigmaPr();

        if (std::abs(nSigmaPi) < sigmaMax) {
          if (std::abs(pion->Rapidity()) > yMax) {
            continue;
          }
          histos.fill(HIST("ptDataTOFPion"), pion->Pt());
          if (!hasTpc)
            histos.fill(HIST("DCAxy_data_pions"), track.dcaXY());
        }
        if (std::abs(nSigmaKa) < sigmaMax) {
          if (std::abs(kaon->Rapidity()) > yMax) {
            continue;
          }
          histos.fill(HIST("ptDataTOFKaon"), kaon->Pt());
          if (!hasTpc)
            histos.fill(HIST("DCAxy_data_kaons"), track.dcaXY());
        }

        if (std::abs(nSigmaPr) < sigmaMax) {
          if (std::abs(proton->Rapidity()) > yMax) {
            continue;
          }
          histos.fill(HIST("ptDataTOFProton"), proton->Pt());
          if (!hasTpc)
            histos.fill(HIST("DCAxy_data_protons"), track.dcaXY());
        }
      }
    }
  }

  PROCESS_SWITCH(PtSpectraInclusiveUpc, processSim, "processSim", false);

  PROCESS_SWITCH(PtSpectraInclusiveUpc, processReco, "processReco", true);

  PROCESS_SWITCH(PtSpectraInclusiveUpc, processData, "processData", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PtSpectraInclusiveUpc>(cfgc)};
}
