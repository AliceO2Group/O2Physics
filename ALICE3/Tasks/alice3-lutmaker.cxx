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

/// \author Nicolo' Jacazio <nicolo.jacazio@cern.ch>, CERN
/// \brief Task to extract LUTs for the fast simulation from full simulation
/// \since 27/04/2021

#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>
#include <SimulationDataFormat/MCUtils.h>

struct Alice3LutMaker {
  static constexpr int nSpecies = 8;
  static constexpr int PDGs[nSpecies] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030};
  static constexpr std::string_view hDirPos[nSpecies] = {"PDG_11/", "PDG_13/", "PDG_211/", "PDG_321/", "PDG_2212/", "PDG_1000010020/", "PDG_1000010030/", "PDG_1000020030/"};
  static constexpr std::string_view hDirNeg[nSpecies] = {"PDG_11/", "PDG_13/", "PDG_211/", "PDG_321/", "PDG_2212/", "PDG_1000010020/", "PDG_1000010030/", "PDG_1000020030/"};
  o2::framework::Configurable<bool> addQA{"addQA", false, "Flag to use add QA plots to show the covariance matrix elements"};
  o2::framework::Configurable<bool> selPrim{"selPrim", false, "If true selects primaries, if not select all particles"};
  o2::framework::Configurable<std::vector<int>> enabledPdgs{"enabledPdgs", std::vector<int>{211, -211}, "List of PDGs enabled"};

  o2::framework::Configurable<int> nchBins{"nchBins", 20, "Number of multiplicity bins"};
  o2::framework::Configurable<float> nchMin{"nchMin", 0.5f, "Lower limit in multiplicity"};
  o2::framework::Configurable<float> nchMax{"nchMax", 3.5f, "Upper limit in multiplicity"};
  o2::framework::Configurable<int> nchLog{"nchLog", 1, "Flag to use a logarithmic multiplicity axis, in this case the Nch limits are the expontents"};

  o2::framework::Configurable<int> etaBins{"etaBins", 80, "Number of eta bins"};
  o2::framework::Configurable<float> etaMin{"etaMin", -4.f, "Lower limit in eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 4.f, "Upper limit in eta"};

  o2::framework::Configurable<int> ptBins{"ptBins", 200, "Number of pT bins"};
  o2::framework::Configurable<float> ptMin{"ptMin", -2.f, "Lower limit in pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 2.f, "Upper limit in pT"};
  o2::framework::Configurable<int> ptLog{"ptLog", 1, "Flag to use a logarithmic pT axis, in this case the pT limits are the expontents"};

  o2::framework::HistogramRegistry histos{"Histos", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};

  void initParticle(o2::track::PID::ID particle, int sign = 1)
  {
    const int pdg = sign * PDGs[particle];
    if (std::find(enabledPdgs.value.begin(), enabledPdgs.value.end(), pdg) == enabledPdgs.value.end()) {
      return;
    }
    LOG(info) << "Initializing LUT maker for PDG " << pdg;
    const TString commonTitle = Form(" PDG %i", pdg);
    const TString pdgdir = (sign > 0 ? hDirPos[particle] : hDirNeg[particle]).data();
    o2::framework::AxisSpec axisPt{ptBins, ptMin, ptMax, "#it{p}_{T} GeV/#it{c}"};
    if (ptLog) {
      if (axisPt.binEdges.size() > 2) {
        LOG(fatal) << "Cannot make a variabled bin width axis logaritmic";
      }
      const double min = axisPt.binEdges[0];
      const double max = axisPt.binEdges[1];
      axisPt.binEdges.clear();
      const int nbins = axisPt.nBins.value();
      const double width = (max - min) / nbins;
      for (int bin = 0; bin < nbins + 1; bin++) {
        float val = min + (bin + 0.5) * width;
        axisPt.binEdges.push_back(pow(10., val));
      }
      axisPt.nBins = std::nullopt;
    }
    o2::framework::AxisSpec axisNch{nchBins, nchMin, nchMax, "N_{Ch}"};
    if (nchLog) {
      axisNch.makeLogarithmic();
    }
    const o2::framework::AxisSpec axisEta{etaBins, etaMin, etaMax, "#it{#eta}"};

    // Track quantities
    histos.add(pdgdir + "multiplicity", "Track multiplicity;Tracks per event;Events", o2::framework::HistType::kTH1F, {axisNch});
    histos.add(pdgdir + "pt", "pt" + commonTitle, o2::framework::HistType::kTH1F, {axisPt});
    histos.add(pdgdir + "eta", "eta" + commonTitle, o2::framework::HistType::kTH1F, {axisEta});

    // Track covariance matrix quantities
    histos.add(pdgdir + "CovMat_sigmaY", "sigmaY" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_sigmaZ", "sigmaZ" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_sigmaSnp", "sigmaSnp" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_sigmaTgl", "sigmaTgl" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_sigma1Pt", "sigma1Pt" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_rhoZY", "rhoZY" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_rhoSnpY", "rhoSnpY" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_rhoSnpZ", "rhoSnpZ" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_rhoTglY", "rhoTglY" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_rhoTglZ", "rhoTglZ" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_rhoTglSnp", "rhoTglSnp" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_rho1PtY", "rho1PtY" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_rho1PtZ", "rho1PtZ" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_rho1PtSnp", "rho1PtSnp" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_rho1PtTgl", "rho1PtTgl" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});

    histos.add(pdgdir + "CovMat_cYY", "cYY" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_cZY", "cZY" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_cZZ", "cZZ" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_cSnpY", "cSnpY" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_cSnpZ", "cSnpZ" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_cSnpSnp", "cSnpSnp" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_cTglY", "cTglY" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_cTglZ", "cTglZ" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_cTglSnp", "cTglSnp" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_cTglTgl", "cTglTgl" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_c1PtY", "c1PtY" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_c1PtZ", "c1PtZ" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_c1PtSnp", "c1PtSnp" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_c1PtTgl", "c1PtTgl" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});
    histos.add(pdgdir + "CovMat_c1Pt21Pt2", "c1Pt21Pt2" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});

    histos.add(pdgdir + "Efficiency", "Efficiency" + commonTitle, o2::framework::HistType::kTProfile2D, {axisPt, axisEta});

    if (!addQA) { // Only if QA histograms are enabled
      return;
    }

    const o2::framework::AxisSpec axissigmaY{300, 0.0, 0.04, "sigmaY"};
    const o2::framework::AxisSpec axissigmaZ{300, 0.0, 0.07, "sigmaZ"};
    const o2::framework::AxisSpec axissigmaSnp{300, 0.0, 0.09, "sigmaSnp"};
    const o2::framework::AxisSpec axissigmaTgl{300, 0.0, 0.15, "sigmaTgl"};
    const o2::framework::AxisSpec axissigma1Pt{300, 0.0, 5, "sigma1Pt"};
    const o2::framework::AxisSpec axisrhoZY{300, -18.0, 18.0, "rhoZY"};
    const o2::framework::AxisSpec axisrhoSnpY{300, -300.0, 0.0, "rhoSnpY"};
    const o2::framework::AxisSpec axisrhoSnpZ{300, -10.0, 10.0, "rhoSnpZ"};
    const o2::framework::AxisSpec axisrhoTglY{300, -20.0, 20.0, "rhoTglY"};
    const o2::framework::AxisSpec axisrhoTglZ{300, -300.0, 0.0, "rhoTglZ"};
    const o2::framework::AxisSpec axisrhoTglSnp{300, -15.0, 15.0, "rhoTglSnp"};
    const o2::framework::AxisSpec axisrho1PtY{321, -250.0, 250.0, "rho1PtY"};
    const o2::framework::AxisSpec axisrho1PtZ{300, -90.0, 90.0, "rho1PtZ"};
    const o2::framework::AxisSpec axisrho1PtSnp{375, -300.0, 300.0, "rho1PtSnp"};
    const o2::framework::AxisSpec axisrho1PtTgl{300, -90.0, 90.0, "rho1PtTgl"};

    histos.add(pdgdir + "QA/CovMat_sigmaY", "sigmaY" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axissigmaY});
    histos.add(pdgdir + "QA/CovMat_sigmaZ", "sigmaZ" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axissigmaZ});
    histos.add(pdgdir + "QA/CovMat_sigmaSnp", "sigmaSnp" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axissigmaSnp});
    histos.add(pdgdir + "QA/CovMat_sigmaTgl", "sigmaTgl" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axissigmaTgl});
    histos.add(pdgdir + "QA/CovMat_sigma1Pt", "sigma1Pt" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axissigma1Pt});
    histos.add(pdgdir + "QA/sigma1Pt", "sigma1Pt" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axissigma1Pt});
    histos.add(pdgdir + "QA/CovMat_rhoZY", "rhoZY" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisrhoZY});
    histos.add(pdgdir + "QA/CovMat_rhoSnpY", "rhoSnpY" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisrhoSnpY});
    histos.add(pdgdir + "QA/CovMat_rhoSnpZ", "rhoSnpZ" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisrhoSnpZ});
    histos.add(pdgdir + "QA/CovMat_rhoTglY", "rhoTglY" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisrhoTglY});
    histos.add(pdgdir + "QA/CovMat_rhoTglZ", "rhoTglZ" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisrhoTglZ});
    histos.add(pdgdir + "QA/CovMat_rhoTglSnp", "rhoTglSnp" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisrhoTglSnp});
    histos.add(pdgdir + "QA/CovMat_rho1PtY", "rho1PtY" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisrho1PtY});
    histos.add(pdgdir + "QA/CovMat_rho1PtZ", "rho1PtZ" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisrho1PtZ});
    histos.add(pdgdir + "QA/CovMat_rho1PtSnp", "rho1PtSnp" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisrho1PtSnp});
    histos.add(pdgdir + "QA/CovMat_rho1PtTgl", "rho1PtTgl" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisrho1PtTgl});

    const o2::framework::AxisSpec axiscYY{300, 0.0, 0.0009, "cYY"};
    const o2::framework::AxisSpec axiscZY{300, -6e-08, 6e-08, "cZY"};
    const o2::framework::AxisSpec axiscZZ{300, 0.0, 0.003, "cZZ"};
    const o2::framework::AxisSpec axiscSnpY{300, -0.0021, 0.0, "cSnpY"};
    const o2::framework::AxisSpec axiscSnpZ{300, -8e-08, 8e-08, "cSnpZ"};
    const o2::framework::AxisSpec axiscSnpSnp{300, 0.0, 0.0025, "cSnpSnp"};
    const o2::framework::AxisSpec axiscTglY{300, -2e-07, 2e-07, "cTglY"};
    const o2::framework::AxisSpec axiscTglZ{300, -0.004, 0.0, "cTglZ"};
    const o2::framework::AxisSpec axiscTglSnp{300, -3.5e-05, 3.5e-05, "cTglSnp"};
    const o2::framework::AxisSpec axiscTglTgl{300, 0.0, 0.008, "cTglTgl"};
    const o2::framework::AxisSpec axisc1PtY{300, -0.004, 0.004, "c1PtY"};
    const o2::framework::AxisSpec axisc1PtZ{300, -0.03, 0.03, "c1PtZ"};
    const o2::framework::AxisSpec axisc1PtSnp{300, -0.015, 0.015, "c1PtSnp"};
    const o2::framework::AxisSpec axisc1PtTgl{300, -0.06, 0.06, "c1PtTgl"};
    const o2::framework::AxisSpec axisc1Pt21Pt2{300, 0.0, 10, "c1Pt21Pt2"};

    histos.add(pdgdir + "QA/CovMat_cYY", "cYY" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axiscYY});
    histos.add(pdgdir + "QA/CovMat_cZY", "cZY" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axiscZY});
    histos.add(pdgdir + "QA/CovMat_cZZ", "cZZ" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axiscZZ});
    histos.add(pdgdir + "QA/CovMat_cSnpY", "cSnpY" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axiscSnpY});
    histos.add(pdgdir + "QA/CovMat_cSnpZ", "cSnpZ" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axiscSnpZ});
    histos.add(pdgdir + "QA/CovMat_cSnpSnp", "cSnpSnp" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axiscSnpSnp});
    histos.add(pdgdir + "QA/CovMat_cTglY", "cTglY" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axiscTglY});
    histos.add(pdgdir + "QA/CovMat_cTglZ", "cTglZ" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axiscTglZ});
    histos.add(pdgdir + "QA/CovMat_cTglSnp", "cTglSnp" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axiscTglSnp});
    histos.add(pdgdir + "QA/CovMat_cTglTgl", "cTglTgl" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axiscTglTgl});
    histos.add(pdgdir + "QA/CovMat_c1PtY", "c1PtY" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisc1PtY});
    histos.add(pdgdir + "QA/CovMat_c1PtZ", "c1PtZ" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisc1PtZ});
    histos.add(pdgdir + "QA/CovMat_c1PtSnp", "c1PtSnp" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisc1PtSnp});
    histos.add(pdgdir + "QA/CovMat_c1PtTgl", "c1PtTgl" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisc1PtTgl});
    histos.add(pdgdir + "QA/CovMat_c1Pt21Pt2", "c1Pt21Pt2" + commonTitle, o2::framework::HistType::kTH3F, {axisPt, axisEta, axisc1Pt21Pt2});
  }

  void init(o2::framework::InitContext&)
  {
    initParticle(o2::track::PID::Pion, 1);
    // initParticle(o2::track::PID::Electron, 1);
    return; // For now only filling the LUT for electrons, to save time, will enable the others later
    for (int i = -1; i >= 1; i += 2) {
      initParticle(o2::track::PID::Electron, i);
      initParticle(o2::track::PID::Muon, i);
      initParticle(o2::track::PID::Kaon, i);
      initParticle(o2::track::PID::Proton, i);
      initParticle(o2::track::PID::Deuteron, i);
      initParticle(o2::track::PID::Triton, i);
      initParticle(o2::track::PID::Helium3, i);
    }
  }

  template <o2::track::PID::ID particle, int sign>
  void processParticle(const o2::aod::McParticles& mcParticles,
                       const o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>::iterator&,
                       const o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::McTrackLabels>& tracks,
                       const o2::aod::McCollisions::iterator&)
  {
    constexpr int pdg = (sign > 0 ? 1 : -1) * PDGs[particle];
    constexpr auto hDIR = HIST((sign > 0 ? hDirPos[particle] : hDirNeg[particle]).data());

    if (std::find(enabledPdgs.value.begin(), enabledPdgs.value.end(), pdg) == enabledPdgs.value.end()) {
      return;
    }
    LOG(info) << "Processing LUT maker for PDG " << pdg;
    std::vector<int64_t> recoTracks(tracks.size());
    int ntrks = 0;

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto& mcParticle = track.mcParticle_as<o2::aod::McParticles>();
      if (mcParticle.pdgCode() != pdg) {
        continue;
      }
      if (selPrim.value && !mcParticle.isPhysicalPrimary()) { // Requiring is physical primary
        continue;
      }

      recoTracks[ntrks++] = mcParticle.globalIndex();

      histos.fill(hDIR + HIST("pt"), mcParticle.pt());
      histos.fill(hDIR + HIST("eta"), mcParticle.eta());

      histos.fill(hDIR + HIST("CovMat_sigmaY"), mcParticle.pt(), mcParticle.eta(), track.sigmaY());
      histos.fill(hDIR + HIST("CovMat_sigmaZ"), mcParticle.pt(), mcParticle.eta(), track.sigmaZ());
      histos.fill(hDIR + HIST("CovMat_sigmaSnp"), mcParticle.pt(), mcParticle.eta(), track.sigmaSnp());
      histos.fill(hDIR + HIST("CovMat_sigmaTgl"), mcParticle.pt(), mcParticle.eta(), track.sigmaTgl());
      histos.fill(hDIR + HIST("CovMat_sigma1Pt"), mcParticle.pt(), mcParticle.eta(), track.sigma1Pt());
      histos.fill(hDIR + HIST("CovMat_rhoZY"), mcParticle.pt(), mcParticle.eta(), track.rhoZY());
      histos.fill(hDIR + HIST("CovMat_rhoSnpY"), mcParticle.pt(), mcParticle.eta(), track.rhoSnpY());
      histos.fill(hDIR + HIST("CovMat_rhoSnpZ"), mcParticle.pt(), mcParticle.eta(), track.rhoSnpZ());
      histos.fill(hDIR + HIST("CovMat_rhoTglY"), mcParticle.pt(), mcParticle.eta(), track.rhoTglY());
      histos.fill(hDIR + HIST("CovMat_rhoTglZ"), mcParticle.pt(), mcParticle.eta(), track.rhoTglZ());
      histos.fill(hDIR + HIST("CovMat_rhoTglSnp"), mcParticle.pt(), mcParticle.eta(), track.rhoTglSnp());
      histos.fill(hDIR + HIST("CovMat_rho1PtY"), mcParticle.pt(), mcParticle.eta(), track.rho1PtY());
      histos.fill(hDIR + HIST("CovMat_rho1PtZ"), mcParticle.pt(), mcParticle.eta(), track.rho1PtZ());
      histos.fill(hDIR + HIST("CovMat_rho1PtSnp"), mcParticle.pt(), mcParticle.eta(), track.rho1PtSnp());
      histos.fill(hDIR + HIST("CovMat_rho1PtTgl"), mcParticle.pt(), mcParticle.eta(), track.rho1PtTgl());

      histos.fill(hDIR + HIST("CovMat_cYY"), mcParticle.pt(), mcParticle.eta(), track.cYY());
      histos.fill(hDIR + HIST("CovMat_cZY"), mcParticle.pt(), mcParticle.eta(), track.cZY());
      histos.fill(hDIR + HIST("CovMat_cZZ"), mcParticle.pt(), mcParticle.eta(), track.cZZ());
      histos.fill(hDIR + HIST("CovMat_cSnpY"), mcParticle.pt(), mcParticle.eta(), track.cSnpY());
      histos.fill(hDIR + HIST("CovMat_cSnpZ"), mcParticle.pt(), mcParticle.eta(), track.cSnpZ());
      histos.fill(hDIR + HIST("CovMat_cSnpSnp"), mcParticle.pt(), mcParticle.eta(), track.cSnpSnp());
      histos.fill(hDIR + HIST("CovMat_cTglY"), mcParticle.pt(), mcParticle.eta(), track.cTglY());
      histos.fill(hDIR + HIST("CovMat_cTglZ"), mcParticle.pt(), mcParticle.eta(), track.cTglZ());
      histos.fill(hDIR + HIST("CovMat_cTglSnp"), mcParticle.pt(), mcParticle.eta(), track.cTglSnp());
      histos.fill(hDIR + HIST("CovMat_cTglTgl"), mcParticle.pt(), mcParticle.eta(), track.cTglTgl());
      histos.fill(hDIR + HIST("CovMat_c1PtY"), mcParticle.pt(), mcParticle.eta(), track.c1PtY());
      histos.fill(hDIR + HIST("CovMat_c1PtZ"), mcParticle.pt(), mcParticle.eta(), track.c1PtZ());
      histos.fill(hDIR + HIST("CovMat_c1PtSnp"), mcParticle.pt(), mcParticle.eta(), track.c1PtSnp());
      histos.fill(hDIR + HIST("CovMat_c1PtTgl"), mcParticle.pt(), mcParticle.eta(), track.c1PtTgl());
      histos.fill(hDIR + HIST("CovMat_c1Pt21Pt2"), mcParticle.pt(), mcParticle.eta(), track.c1Pt21Pt2());

      if (!addQA) { // Only if QA histograms are enabled
        continue;
      }

      histos.fill(hDIR + HIST("QA/CovMat_sigmaY"), mcParticle.pt(), mcParticle.eta(), track.sigmaY());
      histos.fill(hDIR + HIST("QA/CovMat_sigmaZ"), mcParticle.pt(), mcParticle.eta(), track.sigmaZ());
      histos.fill(hDIR + HIST("QA/CovMat_sigmaSnp"), mcParticle.pt(), mcParticle.eta(), track.sigmaSnp());
      histos.fill(hDIR + HIST("QA/CovMat_sigmaTgl"), mcParticle.pt(), mcParticle.eta(), track.sigmaTgl());
      histos.fill(hDIR + HIST("QA/CovMat_sigma1Pt"), mcParticle.pt(), mcParticle.eta(), track.sigma1Pt());
      histos.fill(hDIR + HIST("QA/sigma1Pt"), mcParticle.pt(), mcParticle.eta(), std::abs(track.signed1Pt()) - 1. / mcParticle.pt());
      histos.fill(hDIR + HIST("QA/CovMat_rhoZY"), mcParticle.pt(), mcParticle.eta(), track.rhoZY());
      histos.fill(hDIR + HIST("QA/CovMat_rhoSnpY"), mcParticle.pt(), mcParticle.eta(), track.rhoSnpY());
      histos.fill(hDIR + HIST("QA/CovMat_rhoSnpZ"), mcParticle.pt(), mcParticle.eta(), track.rhoSnpZ());
      histos.fill(hDIR + HIST("QA/CovMat_rhoTglY"), mcParticle.pt(), mcParticle.eta(), track.rhoTglY());
      histos.fill(hDIR + HIST("QA/CovMat_rhoTglZ"), mcParticle.pt(), mcParticle.eta(), track.rhoTglZ());
      histos.fill(hDIR + HIST("QA/CovMat_rhoTglSnp"), mcParticle.pt(), mcParticle.eta(), track.rhoTglSnp());
      histos.fill(hDIR + HIST("QA/CovMat_rho1PtY"), mcParticle.pt(), mcParticle.eta(), track.rho1PtY());
      histos.fill(hDIR + HIST("QA/CovMat_rho1PtZ"), mcParticle.pt(), mcParticle.eta(), track.rho1PtZ());
      histos.fill(hDIR + HIST("QA/CovMat_rho1PtSnp"), mcParticle.pt(), mcParticle.eta(), track.rho1PtSnp());
      histos.fill(hDIR + HIST("QA/CovMat_rho1PtTgl"), mcParticle.pt(), mcParticle.eta(), track.rho1PtTgl());

      histos.fill(hDIR + HIST("QA/CovMat_cYY"), mcParticle.pt(), mcParticle.eta(), track.cYY());
      histos.fill(hDIR + HIST("QA/CovMat_cZY"), mcParticle.pt(), mcParticle.eta(), track.cZY());
      histos.fill(hDIR + HIST("QA/CovMat_cZZ"), mcParticle.pt(), mcParticle.eta(), track.cZZ());
      histos.fill(hDIR + HIST("QA/CovMat_cSnpY"), mcParticle.pt(), mcParticle.eta(), track.cSnpY());
      histos.fill(hDIR + HIST("QA/CovMat_cSnpZ"), mcParticle.pt(), mcParticle.eta(), track.cSnpZ());
      histos.fill(hDIR + HIST("QA/CovMat_cSnpSnp"), mcParticle.pt(), mcParticle.eta(), track.cSnpSnp());
      histos.fill(hDIR + HIST("QA/CovMat_cTglY"), mcParticle.pt(), mcParticle.eta(), track.cTglY());
      histos.fill(hDIR + HIST("QA/CovMat_cTglZ"), mcParticle.pt(), mcParticle.eta(), track.cTglZ());
      histos.fill(hDIR + HIST("QA/CovMat_cTglSnp"), mcParticle.pt(), mcParticle.eta(), track.cTglSnp());
      histos.fill(hDIR + HIST("QA/CovMat_cTglTgl"), mcParticle.pt(), mcParticle.eta(), track.cTglTgl());
      histos.fill(hDIR + HIST("QA/CovMat_c1PtY"), mcParticle.pt(), mcParticle.eta(), track.c1PtY());
      histos.fill(hDIR + HIST("QA/CovMat_c1PtZ"), mcParticle.pt(), mcParticle.eta(), track.c1PtZ());
      histos.fill(hDIR + HIST("QA/CovMat_c1PtSnp"), mcParticle.pt(), mcParticle.eta(), track.c1PtSnp());
      histos.fill(hDIR + HIST("QA/CovMat_c1PtTgl"), mcParticle.pt(), mcParticle.eta(), track.c1PtTgl());
      histos.fill(hDIR + HIST("QA/CovMat_c1Pt21Pt2"), mcParticle.pt(), mcParticle.eta(), track.c1Pt21Pt2());
    }
    histos.fill(hDIR + HIST("multiplicity"), ntrks);

    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.pdgCode() != pdg) {
        continue;
      }
      if (!mcParticle.isPhysicalPrimary()) { // Requiring is physical primary
        continue;
      }

      if (std::find(recoTracks.begin(), recoTracks.end(), mcParticle.globalIndex()) != recoTracks.end()) {
        histos.fill(hDIR + HIST("Efficiency"), mcParticle.pt(), mcParticle.eta(), 1.);
      } else {
        histos.fill(hDIR + HIST("Efficiency"), mcParticle.pt(), mcParticle.eta(), 0.);
      }
    }
  }

  o2::framework::Preslice<o2::aod::McParticles> perMcCollision = o2::aod::mcparticle::mcCollisionId;
  o2::framework::SliceCache cache;
  void process(const o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>::iterator& collision,
               const o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::McTrackLabels>& tracks,
               const o2::aod::McParticles& mcParticles,
               const o2::aod::McCollisions& mcCollisions)
  {
    const auto& mcCollision = collision.mcCollision_as<o2::aod::McCollisions>();
    const auto& particlesInCollision = mcParticles.sliceByCached(o2::aod::mcparticle::mcCollisionId,
                                                                 mcCollision.globalIndex(),
                                                                 cache);
    processParticle<o2::track::PID::Electron, 1>(particlesInCollision, collision, tracks, mcCollision);
    // processParticle<o2::track::PID::Muon, 1>(particlesInCollision, collision, tracks, mcCollision);
    processParticle<o2::track::PID::Pion, 1>(particlesInCollision, collision, tracks, mcCollision);
    // processParticle<o2::track::PID::Kaon, 1>(particlesInCollision, collision, tracks, mcCollision);
    // processParticle<o2::track::PID::Proton, 1>(particlesInCollision, collision, tracks, mcCollision);
    // processParticle<o2::track::PID::Deuteron, 1>(particlesInCollision, collision, tracks, mcCollision);
    // processParticle<o2::track::PID::Triton, 1>(particlesInCollision, collision, tracks, mcCollision);
    // processParticle<o2::track::PID::Helium3, 1>(particlesInCollision, collision, tracks, mcCollision);

    // processParticle<o2::track::PID::Electron, -1>(particlesInCollision, collision, tracks, mcCollision);
    // processParticle<o2::track::PID::Muon, -1>(particlesInCollision, collision, tracks, mcCollision);
    processParticle<o2::track::PID::Pion, -1>(particlesInCollision, collision, tracks, mcCollision);
    // processParticle<o2::track::PID::Kaon, -1>(particlesInCollision, collision, tracks, mcCollision);
    // processParticle<o2::track::PID::Proton, -1>(particlesInCollision, collision, tracks, mcCollision);
    // processParticle<o2::track::PID::Deuteron, -1>(particlesInCollision, collision, tracks, mcCollision);
    // processParticle<o2::track::PID::Triton, -1>(particlesInCollision, collision, tracks, mcCollision);
    // processParticle<o2::track::PID::Helium3, -1>(particlesInCollision, collision, tracks, mcCollision);
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{adaptAnalysisTask<Alice3LutMaker>(cfgc)};
}
