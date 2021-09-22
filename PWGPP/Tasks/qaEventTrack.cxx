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
/// \author Peter Hristov <Peter.Hristov@cern.ch>, CERN
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Henrique J C Zanoli <henrique.zanoli@cern.ch>, Utrecht University
/// \author Nicolo' Jacazio <nicolo.jacazio@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/MC.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2::framework;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{
    {"global", VariantType::Int, 1, {"Run QA of global observable"}},
    {"kine", VariantType::Int, 1, {"Run QA of kinematic observable"}},
    {"reso", VariantType::Int, 1, {"Run QA of resolution observable"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

using namespace o2::framework;
using namespace o2::dataformats;

/// Task to QA global observables of the event
struct QaGlobalObservables {
  // Cuts
  Configurable<int> numberOfContributorsMin{"numberOfContributorsMin", 0, "Minimum required number of contributors to the vertex"};
  Configurable<float> etaMin{"etaMin", -0.8f, "Minimum eta in range to count the track multiplicity"};
  Configurable<float> etaMax{"etaMax", 0.8f, "Maximum eta in range to count the track multiplicity"};

  // Binning
  Configurable<int> numberOfTracksBins{"numberOfTracksBins", 2000, "Number of bins for the Number of Tracks"};
  Configurable<float> numberOfTracksMin{"numberOfTracksMin", 0, "Lower limit in the Number of Tracks plot"};
  Configurable<float> numberOfTracksMax{"numberOfTracksMax", 2000, "Upper limit in the Number of Tracks plot"};

  Configurable<int> vertexPositionBins{"vertexPositionBins", 100, "Number of bins for the Vertex Position"};
  Configurable<float> vertexPositionZMin{"vertexPositionZMin", -20.f, "Lower limit in the Vertex Position Z"};
  Configurable<float> vertexPositionZMax{"vertexPositionZMax", 20.f, "Upper limit in the Vertex Position Z"};
  Configurable<float> vertexPositionXYMin{"vertexPositionXYMin", -0.01f, "Lower limit in the Vertex Position XY"};
  Configurable<float> vertexPositionXYMax{"vertexPositionXYMax", 0.01f, "Upper limit in the Vertex Position XY"};

  Configurable<int> vertexPositionDeltaBins{"vertexPositionDeltaBins", 100, "Number of bins for the histograms of the difference between reconstructed and generated vertex positions"};
  Configurable<float> vertexPositionXYDeltaPtRange{"vertexPositionXYDeltaPtRange", 0.5, "Range of the resolution of the vertex position plot in X and Y"};
  Configurable<float> vertexPositionZDeltaPtRange{"vertexPositionZDeltaPtRange", 0.5, "Range of the resolution of the vertex position plot in Z"};

  Configurable<int> numbersOfContributorsToPVBins{"numbersOfContributorsToPVBins", 200, "Number bins for the number of contributors to the primary vertex"};
  Configurable<float> numbersOfContributorsToPVMax{"numbersOfContributorsToPVMax", 200, "Maximum value for the Number of contributors to the primary vertex"};

  Configurable<int> vertexCovarianceMatrixBins{"vertexCovarianceMatrixBins", 100, "Number bins for the vertex covariance matrix"};
  Configurable<float> vertexCovarianceMatrixMin{"vertexCovarianceMatrixMin", -0.01f, "Lower limit in the Vertex Covariance matrix XY"};
  Configurable<float> vertexCovarianceMatrixMax{"vertexCovarianceMatrixMax", 0.01f, "Upper limit in the Vertex Covariance matrix XY"};

  HistogramRegistry histograms{"HistogramsGlobalQA"};
  void init(InitContext&)
  {

    const AxisSpec numberOfTrackAxis{numberOfTracksBins, numberOfTracksMin, numberOfTracksMax, "Track Multiplicity"};
    const AxisSpec collisionXAxis{vertexPositionBins, vertexPositionXYMin, vertexPositionXYMax, "X [cm]"};
    const AxisSpec collisionYAxis{vertexPositionBins, vertexPositionXYMin, vertexPositionXYMax, "Y [cm]"};
    const AxisSpec collisionZAxis{vertexPositionBins, vertexPositionZMin, vertexPositionZMax, "Z [cm]"};
    const AxisSpec numberOfContributorsAxis{numbersOfContributorsToPVBins, 0, numbersOfContributorsToPVMax, "Number Of contributors to the PV"};
    const AxisSpec vertexCovarianceMatrixAxis{vertexCovarianceMatrixBins, vertexCovarianceMatrixMin, vertexCovarianceMatrixMax};
    const AxisSpec collisionXYDeltaAxis{vertexPositionDeltaBins, -vertexPositionXYDeltaPtRange, vertexPositionXYDeltaPtRange};
    const AxisSpec collisionZDeltaAxis{vertexPositionDeltaBins, -vertexPositionZDeltaPtRange, vertexPositionZDeltaPtRange};

    // Global
    histograms.add("eventCount", ";Selected Events", kTH1D, {{2, 0, 2}});

    // Collision
    histograms.add("collision/Eff", "", kTH1D, {{2, 0.5, 2.5}});
    histograms.get<TH1>(HIST("collision/Eff"))->GetXaxis()->SetBinLabel(1, "Coll. Read");
    histograms.get<TH1>(HIST("collision/Eff"))->GetXaxis()->SetBinLabel(2, "Coll. Reco");
    histograms.add("collision/X", "", kTH1D, {collisionXAxis});
    histograms.add("collision/Y", "", kTH1D, {collisionYAxis});
    histograms.add("collision/Z", "", kTH1D, {collisionZAxis});
    histograms.add("collision/XvsNContrib", "", kTH2D, {collisionXAxis, numberOfContributorsAxis});
    histograms.add("collision/YvsNContrib", "", kTH2D, {collisionYAxis, numberOfContributorsAxis});
    histograms.add("collision/ZvsNContrib", "", kTH2D, {collisionZAxis, numberOfContributorsAxis});
    histograms.add("collision/numberOfContributors", "", kTH1D, {numberOfContributorsAxis});
    histograms.add("collision/numberOfContributorsVsMult", "", kTH2D, {numberOfContributorsAxis, numberOfTrackAxis});
    histograms.add("collision/vertexChi2", ";#chi^{2}", kTH1D, {{100, 0, 10}});

    // Covariance
    histograms.add("covariance/xx", ";Cov_{xx} [cm^{2}]", kTH1D, {vertexCovarianceMatrixAxis});
    histograms.add("covariance/xy", ";Cov_{xy} [cm^{2}]", kTH1D, {vertexCovarianceMatrixAxis});
    histograms.add("covariance/xz", ";Cov_{xz} [cm^{2}]", kTH1D, {vertexCovarianceMatrixAxis});
    histograms.add("covariance/yy", ";Cov_{yy} [cm^{2}]", kTH1D, {vertexCovarianceMatrixAxis});
    histograms.add("covariance/yz", ";Cov_{yz} [cm^{2}]", kTH1D, {vertexCovarianceMatrixAxis});
    histograms.add("covariance/zz", ";Cov_{zz} [cm^{2}]", kTH1D, {vertexCovarianceMatrixAxis});
    // Multiplicity
    histograms.add("multiplicity/numberOfTracks", "", kTH1D, {numberOfTrackAxis});
    // Resolution
    histograms.add("resolution/X", ";X_{Rec} - X_{Gen} [cm]", kTH2D, {collisionXYDeltaAxis, numberOfContributorsAxis});
    histograms.add("resolution/Y", ";Y_{Rec} - Y_{Gen} [cm]", kTH2D, {collisionXYDeltaAxis, numberOfContributorsAxis});
    histograms.add("resolution/Z", ";Z_{Rec} - Z_{Gen} [cm]", kTH2D, {collisionZDeltaAxis, numberOfContributorsAxis});
  }

  void process(const o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>::iterator& collision,
               const o2::aod::McCollisions&,
               const o2::aod::Tracks& tracks)
  {
    histograms.fill(HIST("collision/Eff"), 1);
    if (collision.numContrib() > 0) {
      histograms.fill(HIST("collision/Eff"), 2);
    }

    if (collision.numContrib() < numberOfContributorsMin) {
      return;
    }
    int nTracks = 0;
    for (const auto& track : tracks) {
      if (track.eta() < etaMin || track.eta() > etaMax) {
        continue;
      }
      nTracks++;
    }
    histograms.fill(HIST("eventCount"), 0);

    histograms.fill(HIST("collision/X"), collision.posX());
    histograms.fill(HIST("collision/Y"), collision.posY());
    histograms.fill(HIST("collision/Z"), collision.posZ());

    histograms.fill(HIST("collision/XvsNContrib"), collision.posX(), collision.numContrib());
    histograms.fill(HIST("collision/YvsNContrib"), collision.posY(), collision.numContrib());
    histograms.fill(HIST("collision/ZvsNContrib"), collision.posZ(), collision.numContrib());

    histograms.fill(HIST("collision/numberOfContributors"), collision.numContrib());
    histograms.fill(HIST("collision/numberOfContributorsVsMult"), collision.numContrib(), nTracks);
    histograms.fill(HIST("collision/vertexChi2"), collision.chi2());

    histograms.fill(HIST("covariance/xx"), collision.covXX());
    histograms.fill(HIST("covariance/xy"), collision.covXY());
    histograms.fill(HIST("covariance/xz"), collision.covXZ());
    histograms.fill(HIST("covariance/yy"), collision.covYY());
    histograms.fill(HIST("covariance/yz"), collision.covYZ());
    histograms.fill(HIST("covariance/zz"), collision.covZZ());

    histograms.fill(HIST("multiplicity/numberOfTracks"), nTracks);

    const auto mcColl = collision.mcCollision();
    histograms.fill(HIST("resolution/X"), collision.posX() - mcColl.posX(), collision.numContrib());
    histograms.fill(HIST("resolution/Y"), collision.posY() - mcColl.posY(), collision.numContrib());
    histograms.fill(HIST("resolution/Z"), collision.posZ() - mcColl.posZ(), collision.numContrib());
  }
};

/// Task to QA the kinematic properties of the tracks
struct QaTrackingKine {

  Configurable<int> checkPrimaries{"checkPrimaries", 1,
                                   "Whether to check physical primary and secondaries particles for the resolution."};

  Configurable<int> pdgCodeSel{"pdgCodeSel", 0, "PDG code of the particle to select in absolute value, 0 selects every particle"};

  Configurable<int> ptBins{"ptBins", 100, "Number of pT bins"};
  Configurable<float> ptMin{"ptMin", 0.f, "Lower limit in pT"};
  Configurable<float> ptMax{"ptMax", 5.f, "Upper limit in pT"};

  Configurable<int> etaBins{"etaBins", 100, "Number of eta bins"};
  Configurable<float> etaMin{"etaMin", -3.f, "Lower limit in eta"};
  Configurable<float> etaMax{"etaMax", 3.f, "Upper limit in eta"};

  Configurable<int> phiBins{"phiBins", 100, "Number of phi bins"};
  Configurable<float> phiMin{"phiMin", 0.f, "Lower limit in phi"};
  Configurable<float> phiMax{"phiMax", TMath::TwoPi(), "Upper limit in phi"};

  HistogramRegistry histos{"HistogramsKineQA"};
  void init(InitContext&)
  {
    const AxisSpec ptAxis{ptBins, ptMin, ptMax, "#it{p}_{T} [GeV/#it{c}]"};
    const AxisSpec etaAxis{etaBins, etaMin, etaMax, "#it{#eta}"};
    const AxisSpec phiAxis{phiBins, phiMin, phiMax, "#it{#varphi} [rad]"};

    TString commonTitle = "";
    if (pdgCodeSel != 0) {
      commonTitle += Form("PDG %i", pdgCodeSel.value);
    }

    histos.add("tracking/pt", commonTitle, kTH1D, {ptAxis});
    histos.add("tracking/eta", commonTitle, kTH1D, {etaAxis});
    histos.add("tracking/phi", commonTitle, kTH1D, {phiAxis});

    if (checkPrimaries) {
      histos.add("trackingPrm/pt", commonTitle + " Primary", kTH1D, {ptAxis});
      histos.add("trackingPrm/eta", commonTitle + " Primary", kTH1D, {etaAxis});
      histos.add("trackingPrm/phi", commonTitle + " Primary", kTH1D, {phiAxis});

      histos.add("trackingSec/pt", commonTitle + " Secondary", kTH1D, {ptAxis});
      histos.add("trackingSec/eta", commonTitle + " Secondary", kTH1D, {etaAxis});
      histos.add("trackingSec/phi", commonTitle + " Secondary", kTH1D, {phiAxis});
    }

    histos.add("particle/pt", commonTitle, kTH1D, {ptAxis});
    histos.add("particle/eta", commonTitle, kTH1D, {etaAxis});
    histos.add("particle/phi", commonTitle, kTH1D, {phiAxis});

    if (checkPrimaries) {
      histos.add("particlePrm/pt", commonTitle + " Primary", kTH1D, {ptAxis});
      histos.add("particlePrm/eta", commonTitle + " Primary", kTH1D, {etaAxis});
      histos.add("particlePrm/phi", commonTitle + " Primary", kTH1D, {phiAxis});

      histos.add("particleSec/pt", commonTitle + " Secondary", kTH1D, {ptAxis});
      histos.add("particleSec/eta", commonTitle + " Secondary", kTH1D, {etaAxis});
      histos.add("particleSec/phi", commonTitle + " Secondary", kTH1D, {phiAxis});
    }
  }

  void process(const o2::aod::McParticles& mcParticles,
               const o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::McTrackLabels>& tracks)
  {
    for (const auto& t : tracks) {
      const auto particle = t.mcParticle();
      if (pdgCodeSel != 0 && particle.pdgCode() != pdgCodeSel) { // Checking PDG code
        continue;
      }

      histos.fill(HIST("tracking/pt"), t.pt());
      histos.fill(HIST("tracking/eta"), t.eta());
      histos.fill(HIST("tracking/phi"), t.phi());

      if (!checkPrimaries) {
        continue;
      }

      if (MC::isPhysicalPrimary(particle)) {
        histos.fill(HIST("trackingPrm/pt"), t.pt());
        histos.fill(HIST("trackingPrm/eta"), t.eta());
        histos.fill(HIST("trackingPrm/phi"), t.phi());
      } else {
        histos.fill(HIST("trackingSec/pt"), t.pt());
        histos.fill(HIST("trackingSec/eta"), t.eta());
        histos.fill(HIST("trackingSec/phi"), t.phi());
      }
    }
    for (const auto& particle : mcParticles) {
      if (pdgCodeSel != 0 && particle.pdgCode() != pdgCodeSel) { // Checking PDG code
        continue;
      }
      histos.fill(HIST("particle/pt"), particle.pt());
      histos.fill(HIST("particle/eta"), particle.eta());
      histos.fill(HIST("particle/phi"), particle.phi());

      if (!checkPrimaries) {
        continue;
      }

      if (MC::isPhysicalPrimary(particle)) {
        histos.fill(HIST("particlePrm/pt"), particle.pt());
        histos.fill(HIST("particlePrm/eta"), particle.eta());
        histos.fill(HIST("particlePrm/phi"), particle.phi());
      } else {
        histos.fill(HIST("particleSec/pt"), particle.pt());
        histos.fill(HIST("particleSec/eta"), particle.eta());
        histos.fill(HIST("particleSec/phi"), particle.phi());
      }
    }
  }
};

/// Task to evaluate the tracking resolution (Pt, Eta, Phi and impact parameter)
struct QaTrackingResolution {

  Configurable<int> checkPrimaries{"checkPrimaries", 1,
                                   "Whether to use only physical primary particles for the resolution."};

  Configurable<int> pdgCodeSel{"pdgCodeSel", 0, "PDG code of the particle to select in absolute value, 0 selects every particle"};

  Configurable<int> ptBins{"ptBins", 100, "Number of bins for the transverse momentum"};
  Configurable<float> ptMin{"ptMin", 0.f, "Lower limit in pT"};
  Configurable<float> ptMax{"ptMax", 5.f, "Upper limit in pT"};

  Configurable<int> etaBins{"etaBins", 100, "Number of eta bins"};
  Configurable<float> etaMin{"etaMin", -3.f, "Lower limit in eta"};
  Configurable<float> etaMax{"etaMax", 3.f, "Upper limit in eta"};

  Configurable<int> phiBins{"phiBins", 100, "Number of phi bins"};
  Configurable<float> phiMin{"phiMin", 0.f, "Lower limit in phi"};
  Configurable<float> phiMax{"phiMax", TMath::TwoPi(), "Upper limit in phi"};

  Configurable<int> deltaPtBins{"deltaPtBins", 100, "Number of bins for the transverse momentum differences"};
  Configurable<float> deltaPtMin{"deltaPtMin", -0.5, "Lower limit in delta pT"};
  Configurable<float> deltaPtMax{"deltaPtMax", 0.5, "Upper limit in delta pT"};

  Configurable<int> deltaEtaBins{"deltaEtaBins", 100, "Number of bins for the pseudorapidity differences"};
  Configurable<float> deltaEtaMin{"deltaEtaMin", -0.1, "Lower limit in delta eta"};
  Configurable<float> deltaEtaMax{"deltaEtaMax", 0.1, "Upper limit in delta eta"};

  Configurable<int> deltaPhiBins{"deltaPhiBins", 100, "Number of bins for the azimuthal angle differences"};
  Configurable<float> deltaPhiMin{"deltaPhiMin", -0.1, "Lower limit in delta phi"};
  Configurable<float> deltaPhiMax{"deltaPhiMax", 0.1, "Upper limit in delta phi"};

  Configurable<int> impactParameterBins{"impactParameterBins", 2000, "Number of bins for the Impact parameter"};
  Configurable<float> impactParameterMin{"impactParameterMin", -500, "Lower limit in impact parameter (micrometers)"};
  Configurable<float> impactParameterMax{"impactParameterMax", 500, "Upper limit in impact parameter (micrometers)"};
  Configurable<float> impactParameterResoMin{"impactParameterResoMin", 0, "Lower limit in impact parameter resolution (micrometers)"};
  Configurable<float> impactParameterResoMax{"impactParameterResoMax", 1000, "Upper limit in impact parameter resolution (micrometers)"};

  HistogramRegistry histos{"HistogramsTrackingResolutionQA"};
  void init(InitContext&)
  {
    // Histogram axis definitions

    const AxisSpec ptAxis{ptBins, ptMin, ptMax, "#it{p}_{T} [GeV/#it{c}]"};
    const AxisSpec ptGenAxis{ptBins, ptMin, ptMax, "#it{p}_{T}_{Gen} [GeV/#it{c}]"};
    const AxisSpec ptRecAxis{ptBins, ptMin, ptMax, "#it{p}_{T}_{Rec} [GeV/#it{c}]"};
    const AxisSpec deltaPtAxis{deltaPtBins, deltaPtMin, deltaPtMax, "#it{p}_{T}_{Rec} - #it{p}_{T}_{Gen} [GeV/#it{c}]"};
    const AxisSpec deltaPtRelativeAxis{deltaPtBins, deltaPtMin, deltaPtMax, "(#it{p}_{T}_{Rec} - #it{p}_{T}_{Gen})/(#it{p}_{T}_{Gen})"};

    const AxisSpec etaAxis{etaBins, etaMin, etaMax, "#it{#eta}"};
    const AxisSpec etaGenAxis{etaBins, etaMin, etaMax, "#it{#eta}_{Gen}"};
    const AxisSpec etaRecAxis{etaBins, etaMin, etaMax, "#it{#eta}_{Rec}"};
    const AxisSpec deltaEtaAxis{deltaEtaBins, deltaEtaMin, deltaEtaMax, "#it{#eta}_{Rec} - #it{#eta}_{Gen}"};

    const AxisSpec phiAxis{phiBins, phiMin, phiMax, "#it{#varphi} [rad]"};
    const AxisSpec phiRecAxis{phiBins, phiMin, phiMax, "#it{#varphi}_{Rec} [rad]"};
    const AxisSpec deltaPhiAxis{deltaPhiBins, deltaPhiMin, deltaPhiMax, "#it{#varphi}_{Gen} - #it{#varphi}_{Rec} [rad]"};

    const AxisSpec impactParRPhiAxis{impactParameterBins, impactParameterMin, impactParameterMax, "Impact Parameter r#it{#varphi} [#mum]"};
    const AxisSpec impactParRPhiErrorAxis{impactParameterBins, impactParameterResoMin, impactParameterResoMax, "Impact Parameter Error r#it{#varphi} [#mum]"};

    const AxisSpec impactParZAxis{impactParameterBins, impactParameterMin, impactParameterMax, "Impact Parameter Z [#mum]"};
    const AxisSpec impactParZErrorAxis{impactParameterBins, impactParameterResoMin, impactParameterResoMax, "Impact Parameter Error Z [#mum]"};

    TString commonTitle = "";
    if (pdgCodeSel != 0) {
      commonTitle += Form("PDG %i", pdgCodeSel.value);
    }
    if (checkPrimaries == 1) {
      commonTitle += " Primary";
    }

    // Eta
    histos.add("eta/etaDiffRecGen", commonTitle, kTH1D, {deltaEtaAxis});
    histos.add("eta/etaDiffRecGenVsEtaGen", commonTitle, kTH2D, {deltaEtaAxis, etaGenAxis});
    histos.add("eta/etaDiffRecGenVsEtaRec", commonTitle, kTH2D, {deltaEtaAxis, etaRecAxis});

    // Phi
    histos.add("phi/phiDiffRecGen", commonTitle, kTH1D, {deltaPhiAxis});

    // Pt
    histos.add("pt/ptDiffRecGen", commonTitle, kTH1D, {deltaPtAxis});
    histos.add("pt/ptResolution", commonTitle, kTH1D, {deltaPtRelativeAxis});
    histos.add("pt/ptResolutionVsPt", commonTitle, kTH2D, {ptRecAxis, deltaPtRelativeAxis});
    histos.add("pt/ptResolutionVsEta", commonTitle, kTH2D, {etaAxis, deltaPtRelativeAxis});
    histos.add("pt/ptResolutionVsPhi", commonTitle, kTH2D, {phiAxis, deltaPtRelativeAxis});

    // Impact parameters

    histos.add("impactParameter/impactParameterRPhiVsPt", commonTitle, kTH2D, {ptRecAxis, impactParRPhiAxis});
    histos.add("impactParameter/impactParameterRPhiVsEta", commonTitle, kTH2D, {etaRecAxis, impactParRPhiAxis});
    histos.add("impactParameter/impactParameterRPhiVsPhi", commonTitle, kTH2D, {phiRecAxis, impactParRPhiAxis});

    histos.add("impactParameter/impactParameterErrorRPhiVsPt", commonTitle, kTH2D, {ptRecAxis, impactParRPhiErrorAxis});
    histos.add("impactParameter/impactParameterErrorRPhiVsEta", commonTitle, kTH2D, {etaRecAxis, impactParRPhiErrorAxis});
    histos.add("impactParameter/impactParameterErrorRPhiVsPhi", commonTitle, kTH2D, {phiRecAxis, impactParRPhiErrorAxis});

    histos.add("impactParameter/impactParameterZVsPt", commonTitle, kTH2D, {ptRecAxis, impactParZAxis});
    histos.add("impactParameter/impactParameterZVsEta", commonTitle, kTH2D, {etaRecAxis, impactParZAxis});
    histos.add("impactParameter/impactParameterZVsPhi", commonTitle, kTH2D, {phiRecAxis, impactParZAxis});

    histos.add("impactParameter/impactParameterErrorZVsPt", commonTitle, kTH2D, {ptRecAxis, impactParZErrorAxis});
    histos.add("impactParameter/impactParameterErrorZVsEta", commonTitle, kTH2D, {etaRecAxis, impactParZErrorAxis});
    histos.add("impactParameter/impactParameterErrorZVsPhi", commonTitle, kTH2D, {phiRecAxis, impactParZErrorAxis});
  }

  void process(const o2::aod::McParticles& mcParticles,
               const o2::aod::Collisions& collision,
               const o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::McTrackLabels>& tracks)
  {
    DCA dca;
    // FIXME: get this from CCDB
    constexpr float magneticField{5.0};      // in kG
    constexpr float toMicrometers = 10000.f; // Conversion from [cm] to [mum]
    float impactParameterRPhi = -999.f;
    float impactParameterRPhiError = -999.f;
    float impactParameterZ = -999.f;
    float impactParameterErrorZ = -999.f;

    for (const auto& track : tracks) {

      const auto particle = track.mcParticle();
      if (pdgCodeSel != 0 && particle.pdgCode() != pdgCodeSel) {
        continue;
      }
      if (checkPrimaries && !MC::isPhysicalPrimary(particle)) {
        continue;
      }
      const double deltaPt = track.pt() - particle.pt();
      histos.fill(HIST("pt/ptDiffRecGen"), deltaPt);

      const double deltaPtOverPt = deltaPt / particle.pt();

      histos.fill(HIST("pt/ptResolution"), deltaPtOverPt);
      histos.fill(HIST("pt/ptResolutionVsPt"), track.pt(), deltaPtOverPt);
      histos.fill(HIST("pt/ptResolutionVsEta"), track.eta(), deltaPtOverPt);
      histos.fill(HIST("pt/ptResolutionVsPhi"), track.phi(), deltaPtOverPt);

      const double deltaEta = track.eta() - particle.eta();
      histos.fill(HIST("eta/etaDiffRecGen"), deltaEta);
      histos.fill(HIST("eta/etaDiffRecGenVsEtaGen"), deltaEta, particle.eta());
      histos.fill(HIST("eta/etaDiffRecGenVsEtaRec"), deltaEta, track.eta());

      histos.fill(HIST("phi/phiDiffRecGen"), track.phi() - particle.phi());
      if (getTrackParCov(track).propagateToDCA(getPrimaryVertex(track.collision()), magneticField, &dca, 100.)) { // Check that the propagation is successfull
        impactParameterRPhi = toMicrometers * dca.getY();
        impactParameterRPhiError = toMicrometers * sqrt(dca.getSigmaY2());
        impactParameterZ = toMicrometers * dca.getZ();
        impactParameterErrorZ = toMicrometers * sqrt(dca.getSigmaZ2());

        histos.fill(HIST("impactParameter/impactParameterRPhiVsPt"), track.pt(), impactParameterRPhi);
        histos.fill(HIST("impactParameter/impactParameterRPhiVsEta"), track.eta(), impactParameterRPhi);
        histos.fill(HIST("impactParameter/impactParameterRPhiVsPhi"), track.phi(), impactParameterRPhi);

        histos.fill(HIST("impactParameter/impactParameterZVsPt"), track.pt(), impactParameterZ);
        histos.fill(HIST("impactParameter/impactParameterZVsEta"), track.eta(), impactParameterZ);
        histos.fill(HIST("impactParameter/impactParameterZVsPhi"), track.phi(), impactParameterZ);

        histos.fill(HIST("impactParameter/impactParameterErrorRPhiVsPt"), track.pt(), impactParameterRPhiError);
        histos.fill(HIST("impactParameter/impactParameterErrorRPhiVsEta"), track.eta(), impactParameterRPhiError);
        histos.fill(HIST("impactParameter/impactParameterErrorRPhiVsPhi"), track.phi(), impactParameterRPhiError);

        histos.fill(HIST("impactParameter/impactParameterErrorZVsPt"), track.pt(), impactParameterErrorZ);
        histos.fill(HIST("impactParameter/impactParameterErrorZVsEta"), track.eta(), impactParameterErrorZ);
        histos.fill(HIST("impactParameter/impactParameterErrorZVsPhi"), track.phi(), impactParameterErrorZ);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w;
  if (cfgc.options().get<int>("global")) {
    w.push_back(adaptAnalysisTask<QaGlobalObservables>(cfgc));
  }
  if (cfgc.options().get<int>("kine")) {
    w.push_back(adaptAnalysisTask<QaTrackingKine>(cfgc));
  }
  if (cfgc.options().get<int>("reso")) {
    w.push_back(adaptAnalysisTask<QaTrackingResolution>(cfgc));
  }
  return w;
}
