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
/// \brief Task for pT spectra of pions, kaons and protons in inclusive UPC events.
///        Used to obtain the templates for the DCA_xy fits for the primary fractions.
///        Fractions are obtained separately for TPC and TOF spectra: a particle passing both PID
///        selections will contribute to both histograms.
///
/// \author Andrea Giovanni Riffero andrea.giovanni.riffero@cern.ch

#include "PWGUD/DataModel/UDTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TMCProcess.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct PtSpectraInclusiveUpc {

  Service<o2::framework::O2DatabasePDG> pdg;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ConfigurableAxis ptBinning{"ptBinning",
    {VARIABLE_WIDTH, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                     1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                     2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0},
    "#it{p}_{T} (GeV/#it{c})"};

  ConfigurableAxis dcaXYaxis{"dcaXYaxis", {1000, -0.6, 0.6}, "DCA_{xy} (cm) binning"};
  Configurable<bool> applyKineCutsInGen{"applyKineCutsInGen", false, "Apply kinematic cuts in the generated level"};

  Configurable<double> etaMax{"etaMax", 0.9, "Maximum track pseudorapidity"};
  Configurable<double> yMax{"yMax", 0.9, "Maximum particle rapidity"};
  Configurable<double> ptMin{"ptMin", 0.1, "Minimum track transverse momentum (GeV/c)"};
  Configurable<int> nFindableMin{"nFindableMin", 70, "Minimum number of findable TPC clusters"};
  Configurable<double> sigmaMax{"sigmaMax", 3., "Maximum absolute PID n-sigma"};
  Configurable<double> dcaZlimit{"dcaZlimit", 2., "Maximum absolute DCA in z (cm)"};
  Configurable<double> maxChi2TPC{"maxChi2TPC", 4., "Maximum TPC chi2 per cluster"};
  Configurable<double> maxChi2ITS{"maxChi2ITS", 36., "Maximum ITS chi2 per cluster"};

  // define abbreviations
  using CCs = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>;
  using CC = CCs::iterator;
  using CCMCs = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDMcCollsLabels>;
  using CCMC = CCMCs::iterator;
  using TCs = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using TCMCs = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA, aod::UDMcTrackLabels>;

  void init(InitContext const&)
  {

    // axes
    const AxisSpec axisPt{
      ptBinning,
      "#it{p}_{T} (GeV/#it{c})",
      "axisPt"};

    const AxisSpec axisEventCounter{
      2,
      0.5,
      2.5,
      "Event type"};

    const AxisSpec axisDCAxy{
      dcaXYaxis,
      "DCA_{xy} (cm)"};

    // histograms
    // generated pT spectra
    histos.add("ptGeneratedPion", "Generated pions;#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});
    histos.add("ptGeneratedKaon", "Generated kaons;#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});
    histos.add("ptGeneratedProton", "Generated protons;#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});

    // reconstructed pT spectra with TPC PID
    histos.add("ptReconstructedTPCPion", "Reconstructed pions (TPC);#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});
    histos.add("ptReconstructedTPCKaon", "Reconstructed kaons (TPC);#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});
    histos.add("ptReconstructedTPCProton", "Reconstructed protons (TPC);#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});

    // reconstructed pT spectra with TOF PID
    histos.add("ptReconstructedTOFPion", "Reconstructed TOF pions;#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});
    histos.add("ptReconstructedTOFKaon", "Reconstructed TOF kaons;#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});
    histos.add("ptReconstructedTOFProton", "Reconstructed TOF protons;#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});

    // data pT spectra with TPC PID
    histos.add("ptDataTPCPion", "Data pions (TPC);#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});
    histos.add("ptDataTPCKaon", "Data kaons (TPC);#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});
    histos.add("ptDataTPCProton", "Data protons (TPC);#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});

    // data pT spectra with TOF PID
    histos.add("ptDataTOFPion", "Data pions (TOF);#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});
    histos.add("ptDataTOFKaon", "Data kaons (TOF);#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});
    histos.add("ptDataTOFProton", "Data protons (TOF);#it{p}_{T} (GeV/#it{c});Counts", kTH1F, {axisPt});

    histos.add("myEventCounter", "Event counter;Event type;Counts", kTH1F, {axisEventCounter});

    histos.add(
      "DCAxy_TPC_primary_pions",
      "Primary pions (TPC);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TPC_secondary_pions",
      "Secondary pions (TPC);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TPC_primary_kaons",
      "Primary kaons (TPC);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TPC_secondary_kaons",
      "Secondary kaons (TPC);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TPC_primary_protons",
      "Primary protons (TPC);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TPC_secondary_protons",
      "Secondary protons from decays (TPC);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TPC_material_protons",
      "Secondary protons from material (TPC);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TPC_data_pions",
      "Data pion candidates (TPC);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TPC_data_kaons",
      "Data kaon candidates (TPC);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TPC_data_protons",
      "Data proton candidates (TPC);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TOF_primary_pions",
      "Primary pions (TOF);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TOF_secondary_pions",
      "Secondary pions (TOF);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TOF_primary_kaons",
      "Primary kaons (TOF);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TOF_secondary_kaons",
      "Secondary kaons (TOF);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TOF_primary_protons",
      "Primary protons (TOF);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TOF_secondary_protons",
      "Secondary protons from decays (TOF);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TOF_material_protons",
      "Secondary protons from material (TOF);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TOF_data_pions",
      "Data pion candidates (TOF);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TOF_data_kaons",
      "Data kaon candidates (TOF);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});

    histos.add(
      "DCAxy_TOF_data_protons",
      "Data proton candidates (TOF);#it{p}_{T} (GeV/#it{c});DCA_{xy} (cm)",
      kTH2F,
      {axisPt, axisDCAxy});
  }

  void processSim(aod::UDMcCollision const&, aod::UDMcParticles const& mcParticles)
  {

    std::array<float, 3> trackMomentum;

    for (const auto& mcParticle : mcParticles) {
      if (!mcParticle.isPhysicalPrimary())
        continue;

      trackMomentum[0] = mcParticle.px();
      trackMomentum[1] = mcParticle.py();
      trackMomentum[2] = mcParticle.pz();

      if (applyKineCutsInGen) {
        if (std::fabs(RecoDecay::eta(trackMomentum)) > etaMax)
          continue;

        if (std::fabs(RecoDecay::y(trackMomentum, pdg->Mass(mcParticle.pdgCode()))) > yMax)
          continue;

        if (RecoDecay::pt(trackMomentum) < ptMin)
          continue;
      }

      if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus) {
        histos.fill(HIST("ptGeneratedPion"), RecoDecay::pt(trackMomentum));
      }

      if (std::abs(mcParticle.pdgCode()) == PDG_t::kKPlus) {
        histos.fill(HIST("ptGeneratedKaon"), RecoDecay::pt(trackMomentum));
      }

      if (std::abs(mcParticle.pdgCode()) == PDG_t::kProton) {
        histos.fill(HIST("ptGeneratedProton"), RecoDecay::pt(trackMomentum));
      }

      histos.fill(HIST("myEventCounter"), 1); // gen event
    }
  }

  template <bool isMc, typename TTracks, typename TCollision>
  void selectDataAndFillHistos(TCollision const& /*collision*/, TTracks const& tracks)
  {

    double dcaXyLimit = 0;
    bool passDCAxyCut = true;

    auto nSigmaPi = -999.;
    auto nSigmaKa = -999.;
    auto nSigmaPr = -999.;

    std::array<float, 3> trackMomentum;

    for (const auto& track : tracks) {
      passDCAxyCut = true;

      if (!track.isPVContributor()) {
        continue;
      }

      // here TPC selection on findable number of clusters
      if (track.hasTPC()) {
        if (track.tpcNClsFindable() < nFindableMin) {
          continue;
        }
        if (track.tpcChi2NCl() > maxChi2TPC) {
          continue;
        }
      }

      if(track.itsChi2NCl() > maxChi2ITS) {
        continue;
      }

      if (track.pt() < ptMin) {
        continue;
      }

      if (!(std::abs(track.dcaZ()) < dcaZlimit)) {
        continue;
      }

      dcaXyLimit = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
      if ((std::abs(track.dcaXY()) > dcaXyLimit)) {
        passDCAxyCut = false;
      }

      bool isPrimary = false;
      bool isDecay = false;
      if constexpr (isMc) {
        if (!track.has_udMcParticle()) {
          continue;
        }
        const auto mcParticle = track.udMcParticle();
        isPrimary = mcParticle.isPhysicalPrimary();
        isDecay = mcParticle.getProcess() == kPDecay;
      }

      trackMomentum[0] = track.px();
      trackMomentum[1] = track.py();
      trackMomentum[2] = track.pz();

      // compute rapidity for each ma
      // TPC tracks
      if (track.hasTPC()) {
        nSigmaPi = track.tpcNSigmaPi();
        nSigmaKa = track.tpcNSigmaKa();
        nSigmaPr = track.tpcNSigmaPr();

        if (std::abs(nSigmaPi) < sigmaMax) {
          if (std::abs(RecoDecay::y(trackMomentum, o2::constants::physics::MassPionCharged)) < yMax) {
            if constexpr (isMc) {
              if (isPrimary) {
                if(passDCAxyCut) histos.fill(HIST("ptReconstructedTPCPion"), track.pt());
                histos.fill(HIST("DCAxy_TPC_primary_pions"), track.pt(), track.dcaXY());
              } else {
                histos.fill(HIST("DCAxy_TPC_secondary_pions"), track.pt(), track.dcaXY());
              }
            } else {
              if(passDCAxyCut) histos.fill(HIST("ptDataTPCPion"), track.pt());
              histos.fill(HIST("DCAxy_TPC_data_pions"), track.pt(), track.dcaXY());
            }
          }
        }
        if (std::abs(nSigmaKa) < sigmaMax) {
          if (std::abs(RecoDecay::y(trackMomentum, o2::constants::physics::MassKaonCharged)) < yMax) {
            if constexpr (isMc) {
              if (isPrimary) {
                if(passDCAxyCut) histos.fill(HIST("ptReconstructedTPCKaon"), track.pt());
                histos.fill(HIST("DCAxy_TPC_primary_kaons"), track.pt(), track.dcaXY());
              } else {
                histos.fill(HIST("DCAxy_TPC_secondary_kaons"), track.pt(), track.dcaXY());
              }
            } else {
              if(passDCAxyCut) histos.fill(HIST("ptDataTPCKaon"), track.pt());
              histos.fill(HIST("DCAxy_TPC_data_kaons"), track.pt(), track.dcaXY());
            }
          }
        }

        if (std::abs(nSigmaPr) < sigmaMax) {
          if (std::abs(RecoDecay::y(trackMomentum, o2::constants::physics::MassProton)) < yMax) {
            if constexpr (isMc) {
              if (isPrimary) {
                if(passDCAxyCut) histos.fill(HIST("ptReconstructedTPCProton"), track.pt());
                histos.fill(HIST("DCAxy_TPC_primary_protons"), track.pt(), track.dcaXY());
              } else {
                if (isDecay) {
                  histos.fill(HIST("DCAxy_TPC_secondary_protons"), track.pt(), track.dcaXY());
                } else {
                  histos.fill(HIST("DCAxy_TPC_material_protons"), track.pt(), track.dcaXY());
                }
              }
            } else {
              if(passDCAxyCut) histos.fill(HIST("ptDataTPCProton"), track.pt());
              histos.fill(HIST("DCAxy_TPC_data_protons"), track.pt(), track.dcaXY());
            }
          }
        }
      }

      // TOF tracks
      if (track.hasTOF()) {
        nSigmaPi = track.tofNSigmaPi();
        nSigmaKa = track.tofNSigmaKa();
        nSigmaPr = track.tofNSigmaPr();

        if (std::abs(nSigmaPi) < sigmaMax) {
          if (std::abs(RecoDecay::y(trackMomentum, o2::constants::physics::MassPionCharged)) < yMax) {
            if constexpr (isMc) {
              if (isPrimary) {
                if(passDCAxyCut) histos.fill(HIST("ptReconstructedTOFPion"), track.pt());
                histos.fill(HIST("DCAxy_TOF_primary_pions"), track.pt(), track.dcaXY());
              } else {
                histos.fill(HIST("DCAxy_TOF_secondary_pions"), track.pt(), track.dcaXY());
              }
            } else {
              if(passDCAxyCut) histos.fill(HIST("ptDataTOFPion"), track.pt());
              histos.fill(HIST("DCAxy_TOF_data_pions"), track.pt(), track.dcaXY());
            }
          }
        }

        if (std::abs(nSigmaKa) < sigmaMax) {
          if (std::abs(RecoDecay::y(trackMomentum, o2::constants::physics::MassKaonCharged)) < yMax) {
            if constexpr (isMc) {
              if (isPrimary) {
                if(passDCAxyCut) histos.fill(HIST("ptReconstructedTOFKaon"), track.pt());
                histos.fill(HIST("DCAxy_TOF_primary_kaons"), track.pt(), track.dcaXY());
              } else {
                histos.fill(HIST("DCAxy_TOF_secondary_kaons"), track.pt(), track.dcaXY());
              }
            } else {
              if(passDCAxyCut) histos.fill(HIST("ptDataTOFKaon"), track.pt());
              histos.fill(HIST("DCAxy_TOF_data_kaons"), track.pt(), track.dcaXY());
            }
          }
        }

        if (std::abs(nSigmaPr) < sigmaMax) {
          if (std::abs(RecoDecay::y(trackMomentum, o2::constants::physics::MassProton)) < yMax) {
            if constexpr (isMc) {
              if (isPrimary) {
                if(passDCAxyCut) histos.fill(HIST("ptReconstructedTOFProton"), track.pt());
                histos.fill(HIST("DCAxy_TOF_primary_protons"), track.pt(), track.dcaXY());
              } else {
                if (isDecay) {
                  histos.fill(HIST("DCAxy_TOF_secondary_protons"), track.pt(), track.dcaXY());
                } else {
                  histos.fill(HIST("DCAxy_TOF_material_protons"), track.pt(), track.dcaXY());
                }
              }
            } else {
              if(passDCAxyCut) histos.fill(HIST("ptDataTOFProton"), track.pt());
              histos.fill(HIST("DCAxy_TOF_data_protons"), track.pt(), track.dcaXY());
            }
          }
        }

      }
    }
  }

  void processReco(CCMC const& collision, TCMCs const& tracks, aod::UDMcParticles const&)
  {
    selectDataAndFillHistos<true>(collision, tracks);
  }

  void processData(CC const& collision, TCs const& tracks)
  {
    selectDataAndFillHistos<false>(collision, tracks);
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
