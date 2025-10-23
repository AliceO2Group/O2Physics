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

// O2 includes
#include "Framework/AnalysisTask.h"
#include "ReconstructionDataFormats/Track.h"
#include "SimulationDataFormat/MCUtils.h"

using namespace o2;
using namespace framework;
using namespace framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{
    {"lut-el", VariantType::Int, 0, {"LUT input for the Electron PDG code"}},
    {"lut-mu", VariantType::Int, 0, {"LUT input for the Muon PDG code"}},
    {"lut-pi", VariantType::Int, 1, {"LUT input for the Pion PDG code"}},
    {"lut-ka", VariantType::Int, 0, {"LUT input for the Kaon PDG code"}},
    {"lut-pr", VariantType::Int, 0, {"LUT input for the Proton PDG code"}},
    {"lut-tr", VariantType::Int, 0, {"LUT input for the Triton PDG code"}},
    {"lut-de", VariantType::Int, 0, {"LUT input for the Deuteron PDG code"}},
    {"lut-he", VariantType::Int, 0, {"LUT input for the Helium3 PDG code"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

template <o2::track::pid_constants::ID particle>
struct Alice3LutMaker {
  static constexpr int nSpecies = 8;
  static constexpr int PDGs[nSpecies] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030};
  static_assert(particle < nSpecies && "Maximum of particles reached");
  static constexpr int pdg = PDGs[particle];
  Configurable<bool> addQA{"addQA", false, "Flag to use add QA plots to show the covariance matrix elements"};
  Configurable<bool> selPrim{"selPrim", false, "If true selects primaries, if not select all particles"};

  Configurable<int> nchBins{"nchBins", 20, "Number of multiplicity bins"};
  Configurable<float> nchMin{"nchMin", 0.5f, "Lower limit in multiplicity"};
  Configurable<float> nchMax{"nchMax", 3.5f, "Upper limit in multiplicity"};
  Configurable<int> nchLog{"nchLog", 1, "Flag to use a logarithmic multiplicity axis, in this case the Nch limits are the expontents"};

  Configurable<int> etaBins{"etaBins", 80, "Number of eta bins"};
  Configurable<float> etaMin{"etaMin", -4.f, "Lower limit in eta"};
  Configurable<float> etaMax{"etaMax", 4.f, "Upper limit in eta"};

  Configurable<int> ptBins{"ptBins", 200, "Number of pT bins"};
  Configurable<float> ptMin{"ptMin", -2.f, "Lower limit in pT"};
  Configurable<float> ptMax{"ptMax", 2.f, "Upper limit in pT"};
  Configurable<int> ptLog{"ptLog", 1, "Flag to use a logarithmic pT axis, in this case the pT limits are the expontents"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    const TString commonTitle = Form(" PDG %i", pdg);
    AxisSpec axisPt{ptBins, ptMin, ptMax, "#it{p}_{T} GeV/#it{c}"};
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
    AxisSpec axisNch{nchBins, nchMin, nchMax, "N_{Ch}"};
    if (nchLog) {
      axisNch.makeLogarithmic();
    }
    const AxisSpec axisEta{etaBins, etaMin, etaMax, "#it{#eta}"};

    // Track quantities
    histos.add("multiplicity", "Track multiplicity;Tracks per event;Events", kTH1F, {axisNch});
    histos.add("pt", "pt" + commonTitle, kTH1F, {axisPt});
    histos.add("eta", "eta" + commonTitle, kTH1F, {axisEta});

    // Track covariance matrix quantities
    histos.add("CovMat_sigmaY", "sigmaY" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_sigmaZ", "sigmaZ" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_sigmaSnp", "sigmaSnp" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_sigmaTgl", "sigmaTgl" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_sigma1Pt", "sigma1Pt" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_rhoZY", "rhoZY" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_rhoSnpY", "rhoSnpY" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_rhoSnpZ", "rhoSnpZ" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_rhoTglY", "rhoTglY" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_rhoTglZ", "rhoTglZ" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_rhoTglSnp", "rhoTglSnp" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_rho1PtY", "rho1PtY" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_rho1PtZ", "rho1PtZ" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_rho1PtSnp", "rho1PtSnp" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_rho1PtTgl", "rho1PtTgl" + commonTitle, kTProfile2D, {axisPt, axisEta});

    histos.add("CovMat_cYY", "cYY" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_cZY", "cZY" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_cZZ", "cZZ" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_cSnpY", "cSnpY" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_cSnpZ", "cSnpZ" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_cSnpSnp", "cSnpSnp" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_cTglY", "cTglY" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_cTglZ", "cTglZ" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_cTglSnp", "cTglSnp" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_cTglTgl", "cTglTgl" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_c1PtY", "c1PtY" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_c1PtZ", "c1PtZ" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_c1PtSnp", "c1PtSnp" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_c1PtTgl", "c1PtTgl" + commonTitle, kTProfile2D, {axisPt, axisEta});
    histos.add("CovMat_c1Pt21Pt2", "c1Pt21Pt2" + commonTitle, kTProfile2D, {axisPt, axisEta});

    histos.add("Efficiency", "Efficiency" + commonTitle, kTProfile2D, {axisPt, axisEta});

    if (!addQA) { // Only if QA histograms are enabled
      return;
    }

    const AxisSpec axissigmaY{300, 0.0, 0.04, "sigmaY"};
    const AxisSpec axissigmaZ{300, 0.0, 0.07, "sigmaZ"};
    const AxisSpec axissigmaSnp{300, 0.0, 0.09, "sigmaSnp"};
    const AxisSpec axissigmaTgl{300, 0.0, 0.15, "sigmaTgl"};
    const AxisSpec axissigma1Pt{300, 0.0, 5, "sigma1Pt"};
    const AxisSpec axisrhoZY{300, -18.0, 18.0, "rhoZY"};
    const AxisSpec axisrhoSnpY{300, -300.0, 0.0, "rhoSnpY"};
    const AxisSpec axisrhoSnpZ{300, -10.0, 10.0, "rhoSnpZ"};
    const AxisSpec axisrhoTglY{300, -20.0, 20.0, "rhoTglY"};
    const AxisSpec axisrhoTglZ{300, -300.0, 0.0, "rhoTglZ"};
    const AxisSpec axisrhoTglSnp{300, -15.0, 15.0, "rhoTglSnp"};
    const AxisSpec axisrho1PtY{321, -250.0, 250.0, "rho1PtY"};
    const AxisSpec axisrho1PtZ{300, -90.0, 90.0, "rho1PtZ"};
    const AxisSpec axisrho1PtSnp{375, -300.0, 300.0, "rho1PtSnp"};
    const AxisSpec axisrho1PtTgl{300, -90.0, 90.0, "rho1PtTgl"};

    histos.add("QA/CovMat_sigmaY", "sigmaY" + commonTitle, kTH3F, {axisPt, axisEta, axissigmaY});
    histos.add("QA/CovMat_sigmaZ", "sigmaZ" + commonTitle, kTH3F, {axisPt, axisEta, axissigmaZ});
    histos.add("QA/CovMat_sigmaSnp", "sigmaSnp" + commonTitle, kTH3F, {axisPt, axisEta, axissigmaSnp});
    histos.add("QA/CovMat_sigmaTgl", "sigmaTgl" + commonTitle, kTH3F, {axisPt, axisEta, axissigmaTgl});
    histos.add("QA/CovMat_sigma1Pt", "sigma1Pt" + commonTitle, kTH3F, {axisPt, axisEta, axissigma1Pt});
    histos.add("QA/sigma1Pt", "sigma1Pt" + commonTitle, kTH3F, {axisPt, axisEta, axissigma1Pt});
    histos.add("QA/CovMat_rhoZY", "rhoZY" + commonTitle, kTH3F, {axisPt, axisEta, axisrhoZY});
    histos.add("QA/CovMat_rhoSnpY", "rhoSnpY" + commonTitle, kTH3F, {axisPt, axisEta, axisrhoSnpY});
    histos.add("QA/CovMat_rhoSnpZ", "rhoSnpZ" + commonTitle, kTH3F, {axisPt, axisEta, axisrhoSnpZ});
    histos.add("QA/CovMat_rhoTglY", "rhoTglY" + commonTitle, kTH3F, {axisPt, axisEta, axisrhoTglY});
    histos.add("QA/CovMat_rhoTglZ", "rhoTglZ" + commonTitle, kTH3F, {axisPt, axisEta, axisrhoTglZ});
    histos.add("QA/CovMat_rhoTglSnp", "rhoTglSnp" + commonTitle, kTH3F, {axisPt, axisEta, axisrhoTglSnp});
    histos.add("QA/CovMat_rho1PtY", "rho1PtY" + commonTitle, kTH3F, {axisPt, axisEta, axisrho1PtY});
    histos.add("QA/CovMat_rho1PtZ", "rho1PtZ" + commonTitle, kTH3F, {axisPt, axisEta, axisrho1PtZ});
    histos.add("QA/CovMat_rho1PtSnp", "rho1PtSnp" + commonTitle, kTH3F, {axisPt, axisEta, axisrho1PtSnp});
    histos.add("QA/CovMat_rho1PtTgl", "rho1PtTgl" + commonTitle, kTH3F, {axisPt, axisEta, axisrho1PtTgl});

    const AxisSpec axiscYY{300, 0.0, 0.0009, "cYY"};
    const AxisSpec axiscZY{300, -6e-08, 6e-08, "cZY"};
    const AxisSpec axiscZZ{300, 0.0, 0.003, "cZZ"};
    const AxisSpec axiscSnpY{300, -0.0021, 0.0, "cSnpY"};
    const AxisSpec axiscSnpZ{300, -8e-08, 8e-08, "cSnpZ"};
    const AxisSpec axiscSnpSnp{300, 0.0, 0.0025, "cSnpSnp"};
    const AxisSpec axiscTglY{300, -2e-07, 2e-07, "cTglY"};
    const AxisSpec axiscTglZ{300, -0.004, 0.0, "cTglZ"};
    const AxisSpec axiscTglSnp{300, -3.5e-05, 3.5e-05, "cTglSnp"};
    const AxisSpec axiscTglTgl{300, 0.0, 0.008, "cTglTgl"};
    const AxisSpec axisc1PtY{300, -0.004, 0.004, "c1PtY"};
    const AxisSpec axisc1PtZ{300, -0.03, 0.03, "c1PtZ"};
    const AxisSpec axisc1PtSnp{300, -0.015, 0.015, "c1PtSnp"};
    const AxisSpec axisc1PtTgl{300, -0.06, 0.06, "c1PtTgl"};
    const AxisSpec axisc1Pt21Pt2{300, 0.0, 10, "c1Pt21Pt2"};

    histos.add("QA/CovMat_cYY", "cYY" + commonTitle, kTH3F, {axisPt, axisEta, axiscYY});
    histos.add("QA/CovMat_cZY", "cZY" + commonTitle, kTH3F, {axisPt, axisEta, axiscZY});
    histos.add("QA/CovMat_cZZ", "cZZ" + commonTitle, kTH3F, {axisPt, axisEta, axiscZZ});
    histos.add("QA/CovMat_cSnpY", "cSnpY" + commonTitle, kTH3F, {axisPt, axisEta, axiscSnpY});
    histos.add("QA/CovMat_cSnpZ", "cSnpZ" + commonTitle, kTH3F, {axisPt, axisEta, axiscSnpZ});
    histos.add("QA/CovMat_cSnpSnp", "cSnpSnp" + commonTitle, kTH3F, {axisPt, axisEta, axiscSnpSnp});
    histos.add("QA/CovMat_cTglY", "cTglY" + commonTitle, kTH3F, {axisPt, axisEta, axiscTglY});
    histos.add("QA/CovMat_cTglZ", "cTglZ" + commonTitle, kTH3F, {axisPt, axisEta, axiscTglZ});
    histos.add("QA/CovMat_cTglSnp", "cTglSnp" + commonTitle, kTH3F, {axisPt, axisEta, axiscTglSnp});
    histos.add("QA/CovMat_cTglTgl", "cTglTgl" + commonTitle, kTH3F, {axisPt, axisEta, axiscTglTgl});
    histos.add("QA/CovMat_c1PtY", "c1PtY" + commonTitle, kTH3F, {axisPt, axisEta, axisc1PtY});
    histos.add("QA/CovMat_c1PtZ", "c1PtZ" + commonTitle, kTH3F, {axisPt, axisEta, axisc1PtZ});
    histos.add("QA/CovMat_c1PtSnp", "c1PtSnp" + commonTitle, kTH3F, {axisPt, axisEta, axisc1PtSnp});
    histos.add("QA/CovMat_c1PtTgl", "c1PtTgl" + commonTitle, kTH3F, {axisPt, axisEta, axisc1PtTgl});
    histos.add("QA/CovMat_c1Pt21Pt2", "c1Pt21Pt2" + commonTitle, kTH3F, {axisPt, axisEta, axisc1Pt21Pt2});
  }

  void process(const o2::aod::McParticles& mcParticles,
               const o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>&,
               const o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::McTrackLabels>& tracks,
               const o2::aod::McCollisions&)
  {
    std::vector<int64_t> recoTracks(tracks.size());
    int ntrks = 0;

    for (const auto& track : tracks) {
      if (!track.has_mcParticle())
        continue;
      const auto mcParticle = track.mcParticle_as<aod::McParticles>();
      if (mcParticle.pdgCode() != pdg) {
        continue;
      }
      if (selPrim.value && !mcParticle.isPhysicalPrimary()) { // Requiring is physical primary
        continue;
      }

      recoTracks[ntrks++] = mcParticle.globalIndex();

      histos.fill(HIST("pt"), mcParticle.pt());
      histos.fill(HIST("eta"), mcParticle.eta());

      histos.fill(HIST("CovMat_sigmaY"), mcParticle.pt(), mcParticle.eta(), track.sigmaY());
      histos.fill(HIST("CovMat_sigmaZ"), mcParticle.pt(), mcParticle.eta(), track.sigmaZ());
      histos.fill(HIST("CovMat_sigmaSnp"), mcParticle.pt(), mcParticle.eta(), track.sigmaSnp());
      histos.fill(HIST("CovMat_sigmaTgl"), mcParticle.pt(), mcParticle.eta(), track.sigmaTgl());
      histos.fill(HIST("CovMat_sigma1Pt"), mcParticle.pt(), mcParticle.eta(), track.sigma1Pt());
      histos.fill(HIST("CovMat_rhoZY"), mcParticle.pt(), mcParticle.eta(), track.rhoZY());
      histos.fill(HIST("CovMat_rhoSnpY"), mcParticle.pt(), mcParticle.eta(), track.rhoSnpY());
      histos.fill(HIST("CovMat_rhoSnpZ"), mcParticle.pt(), mcParticle.eta(), track.rhoSnpZ());
      histos.fill(HIST("CovMat_rhoTglY"), mcParticle.pt(), mcParticle.eta(), track.rhoTglY());
      histos.fill(HIST("CovMat_rhoTglZ"), mcParticle.pt(), mcParticle.eta(), track.rhoTglZ());
      histos.fill(HIST("CovMat_rhoTglSnp"), mcParticle.pt(), mcParticle.eta(), track.rhoTglSnp());
      histos.fill(HIST("CovMat_rho1PtY"), mcParticle.pt(), mcParticle.eta(), track.rho1PtY());
      histos.fill(HIST("CovMat_rho1PtZ"), mcParticle.pt(), mcParticle.eta(), track.rho1PtZ());
      histos.fill(HIST("CovMat_rho1PtSnp"), mcParticle.pt(), mcParticle.eta(), track.rho1PtSnp());
      histos.fill(HIST("CovMat_rho1PtTgl"), mcParticle.pt(), mcParticle.eta(), track.rho1PtTgl());

      histos.fill(HIST("CovMat_cYY"), mcParticle.pt(), mcParticle.eta(), track.cYY());
      histos.fill(HIST("CovMat_cZY"), mcParticle.pt(), mcParticle.eta(), track.cZY());
      histos.fill(HIST("CovMat_cZZ"), mcParticle.pt(), mcParticle.eta(), track.cZZ());
      histos.fill(HIST("CovMat_cSnpY"), mcParticle.pt(), mcParticle.eta(), track.cSnpY());
      histos.fill(HIST("CovMat_cSnpZ"), mcParticle.pt(), mcParticle.eta(), track.cSnpZ());
      histos.fill(HIST("CovMat_cSnpSnp"), mcParticle.pt(), mcParticle.eta(), track.cSnpSnp());
      histos.fill(HIST("CovMat_cTglY"), mcParticle.pt(), mcParticle.eta(), track.cTglY());
      histos.fill(HIST("CovMat_cTglZ"), mcParticle.pt(), mcParticle.eta(), track.cTglZ());
      histos.fill(HIST("CovMat_cTglSnp"), mcParticle.pt(), mcParticle.eta(), track.cTglSnp());
      histos.fill(HIST("CovMat_cTglTgl"), mcParticle.pt(), mcParticle.eta(), track.cTglTgl());
      histos.fill(HIST("CovMat_c1PtY"), mcParticle.pt(), mcParticle.eta(), track.c1PtY());
      histos.fill(HIST("CovMat_c1PtZ"), mcParticle.pt(), mcParticle.eta(), track.c1PtZ());
      histos.fill(HIST("CovMat_c1PtSnp"), mcParticle.pt(), mcParticle.eta(), track.c1PtSnp());
      histos.fill(HIST("CovMat_c1PtTgl"), mcParticle.pt(), mcParticle.eta(), track.c1PtTgl());
      histos.fill(HIST("CovMat_c1Pt21Pt2"), mcParticle.pt(), mcParticle.eta(), track.c1Pt21Pt2());

      if (!addQA) { // Only if QA histograms are enabled
        continue;
      }

      histos.fill(HIST("QA/CovMat_sigmaY"), mcParticle.pt(), mcParticle.eta(), track.sigmaY());
      histos.fill(HIST("QA/CovMat_sigmaZ"), mcParticle.pt(), mcParticle.eta(), track.sigmaZ());
      histos.fill(HIST("QA/CovMat_sigmaSnp"), mcParticle.pt(), mcParticle.eta(), track.sigmaSnp());
      histos.fill(HIST("QA/CovMat_sigmaTgl"), mcParticle.pt(), mcParticle.eta(), track.sigmaTgl());
      histos.fill(HIST("QA/CovMat_sigma1Pt"), mcParticle.pt(), mcParticle.eta(), track.sigma1Pt());
      histos.fill(HIST("QA/sigma1Pt"), mcParticle.pt(), mcParticle.eta(), std::abs(track.signed1Pt()) - 1. / mcParticle.pt());
      histos.fill(HIST("QA/CovMat_rhoZY"), mcParticle.pt(), mcParticle.eta(), track.rhoZY());
      histos.fill(HIST("QA/CovMat_rhoSnpY"), mcParticle.pt(), mcParticle.eta(), track.rhoSnpY());
      histos.fill(HIST("QA/CovMat_rhoSnpZ"), mcParticle.pt(), mcParticle.eta(), track.rhoSnpZ());
      histos.fill(HIST("QA/CovMat_rhoTglY"), mcParticle.pt(), mcParticle.eta(), track.rhoTglY());
      histos.fill(HIST("QA/CovMat_rhoTglZ"), mcParticle.pt(), mcParticle.eta(), track.rhoTglZ());
      histos.fill(HIST("QA/CovMat_rhoTglSnp"), mcParticle.pt(), mcParticle.eta(), track.rhoTglSnp());
      histos.fill(HIST("QA/CovMat_rho1PtY"), mcParticle.pt(), mcParticle.eta(), track.rho1PtY());
      histos.fill(HIST("QA/CovMat_rho1PtZ"), mcParticle.pt(), mcParticle.eta(), track.rho1PtZ());
      histos.fill(HIST("QA/CovMat_rho1PtSnp"), mcParticle.pt(), mcParticle.eta(), track.rho1PtSnp());
      histos.fill(HIST("QA/CovMat_rho1PtTgl"), mcParticle.pt(), mcParticle.eta(), track.rho1PtTgl());

      histos.fill(HIST("QA/CovMat_cYY"), mcParticle.pt(), mcParticle.eta(), track.cYY());
      histos.fill(HIST("QA/CovMat_cZY"), mcParticle.pt(), mcParticle.eta(), track.cZY());
      histos.fill(HIST("QA/CovMat_cZZ"), mcParticle.pt(), mcParticle.eta(), track.cZZ());
      histos.fill(HIST("QA/CovMat_cSnpY"), mcParticle.pt(), mcParticle.eta(), track.cSnpY());
      histos.fill(HIST("QA/CovMat_cSnpZ"), mcParticle.pt(), mcParticle.eta(), track.cSnpZ());
      histos.fill(HIST("QA/CovMat_cSnpSnp"), mcParticle.pt(), mcParticle.eta(), track.cSnpSnp());
      histos.fill(HIST("QA/CovMat_cTglY"), mcParticle.pt(), mcParticle.eta(), track.cTglY());
      histos.fill(HIST("QA/CovMat_cTglZ"), mcParticle.pt(), mcParticle.eta(), track.cTglZ());
      histos.fill(HIST("QA/CovMat_cTglSnp"), mcParticle.pt(), mcParticle.eta(), track.cTglSnp());
      histos.fill(HIST("QA/CovMat_cTglTgl"), mcParticle.pt(), mcParticle.eta(), track.cTglTgl());
      histos.fill(HIST("QA/CovMat_c1PtY"), mcParticle.pt(), mcParticle.eta(), track.c1PtY());
      histos.fill(HIST("QA/CovMat_c1PtZ"), mcParticle.pt(), mcParticle.eta(), track.c1PtZ());
      histos.fill(HIST("QA/CovMat_c1PtSnp"), mcParticle.pt(), mcParticle.eta(), track.c1PtSnp());
      histos.fill(HIST("QA/CovMat_c1PtTgl"), mcParticle.pt(), mcParticle.eta(), track.c1PtTgl());
      histos.fill(HIST("QA/CovMat_c1Pt21Pt2"), mcParticle.pt(), mcParticle.eta(), track.c1Pt21Pt2());
    }
    histos.fill(HIST("multiplicity"), ntrks);

    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.pdgCode() != pdg) {
        continue;
      }
      if (!mcParticle.isPhysicalPrimary()) { // Requiring is physical primary
        continue;
      }

      if (std::find(recoTracks.begin(), recoTracks.end(), mcParticle.globalIndex()) != recoTracks.end()) {
        histos.fill(HIST("Efficiency"), mcParticle.pt(), mcParticle.eta(), 1.);
      } else {
        histos.fill(HIST("Efficiency"), mcParticle.pt(), mcParticle.eta(), 0.);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w;
  if (cfgc.options().get<int>("lut-el")) {
    w.push_back(adaptAnalysisTask<Alice3LutMaker<o2::track::PID::Electron>>(cfgc, TaskName{"alice3-lutmaker-electron"}));
  }
  if (cfgc.options().get<int>("lut-mu")) {
    w.push_back(adaptAnalysisTask<Alice3LutMaker<o2::track::PID::Muon>>(cfgc, TaskName{"alice3-lutmaker-muon"}));
  }
  if (cfgc.options().get<int>("lut-pi")) {
    w.push_back(adaptAnalysisTask<Alice3LutMaker<o2::track::PID::Pion>>(cfgc, TaskName{"alice3-lutmaker-pion"}));
  }
  if (cfgc.options().get<int>("lut-ka")) {
    w.push_back(adaptAnalysisTask<Alice3LutMaker<o2::track::PID::Kaon>>(cfgc, TaskName{"alice3-lutmaker-kaon"}));
  }
  if (cfgc.options().get<int>("lut-pr")) {
    w.push_back(adaptAnalysisTask<Alice3LutMaker<o2::track::PID::Proton>>(cfgc, TaskName{"alice3-lutmaker-proton"}));
  }
  if (cfgc.options().get<int>("lut-de")) {
    w.push_back(adaptAnalysisTask<Alice3LutMaker<o2::track::PID::Deuteron>>(cfgc, TaskName{"alice3-lutmaker-deuteron"}));
  }
  if (cfgc.options().get<int>("lut-tr")) {
    w.push_back(adaptAnalysisTask<Alice3LutMaker<o2::track::PID::Triton>>(cfgc, TaskName{"alice3-lutmaker-triton"}));
  }
  if (cfgc.options().get<int>("lut-he")) {
    w.push_back(adaptAnalysisTask<Alice3LutMaker<o2::track::PID::Helium3>>(cfgc, TaskName{"alice3-lutmaker-helium3"}));
  }
  return w;
}
