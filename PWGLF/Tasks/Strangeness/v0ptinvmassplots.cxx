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
/// \file v0ptinvmassplots.cxx
/// \brief V0 task for production of invariant mass plots for Pt Spectrum Analysis
/// \author Nikolaos Karatzenis (nikolaos.karatzenis@cern.ch)
/// \author Roman Lietava (roman.lietava@cern.ch)

/*Description
This task creates up to 30 histograms that are filled with the V0 invariant mass under the K0, Lambda and Antilambda mass assumption
for different pt ranges (constituting bins). The values are inserted as configurable strings for convinience.
Also feed-down matrices for the Lambda and Anti-Lambda are produced.
This analysis includes three processes, one for Real Data and two for MC at the Generated and Reconstructed level*/

#include <memory>
#include <vector>
#include <string>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonUtils/StringUtils.h"

// namespace to be used for pt plots and bins
namespace pthistos
{
std::vector<std::shared_ptr<TH1>> kaonPt;
static std::vector<std::string> kaonPtBins;
std::vector<std::shared_ptr<TH1>> lambdaPt;
static std::vector<std::string> lambdaPtBins;
std::vector<std::shared_ptr<TH1>> antilambdaPt;
static std::vector<std::string> antilambdaPtBins;
} // namespace pthistos
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct V0PtInvMassPlots {
  // Histogram Registries
  HistogramRegistry rPtAnalysis{"PtAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKaonshMassPlotsPerPtBin{"KaonshMassPlotsPerPtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambdaMassPlotsPerPtBin{"LambdaMassPlotsPerPtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntilambdaMassPlotsPerPtBin{"AntilambdaMassPlotsPerPtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rFeeddownMatrices{"FeeddownMatrices", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  Configurable<int> nBinsArmenteros{"nBinsArmenteros", 500, "N bins in Armenteros histos"};
  Configurable<int> nmaxHistograms{"nmaxHistograms", 20, "N Pt Histograms"};

  // Configurables for Cuts
  Configurable<float> cutZVertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> nSigmaTPCPion{"nSigmaTPCPion", 4, "nSigmaTPCPion"};
  Configurable<float> nSigmaTPCProton{"nSigmaTPCProton", 4, "nSigmaTPCProton"};
  Configurable<float> compv0masscut{"compv0masscut", 0.01, "CompetitiveV0masscut (GeV)"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};
  Configurable<float> rapiditycutGen{"rapiditycutGen", 0.5, "V0 Rapidity Window GenMC"};

  // Configurable Kaonsh Topological Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> kaonshSettingdcav0dau{"kaonshSettingdcav0dau", 0.3, "DCA V0 Daughters"};
  Configurable<float> kaonshSettingdcapostopv{"kaonshSettingdcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<float> kaonshSettingdcanegtopv{"kaonshSettingdcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<double> kaonshSettingcosPA{"kaonshSettingcosPA", 0.98, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> kaonshSettingradius{"kaonshSettingradius", 0.50, "v0radius"};

  // Configurable Lambda Topological Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> lambdaSettingdcav0dau{"lambdaSettingdcav0dau", 0.3, "DCA V0 Daughters"};
  Configurable<float> lambdaSettingdcapostopv{"lambdaSettingdcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<float> lambdaSettingdcanegtopv{"lambdaSettingdcanegtopv", 0.09, "DCA Neg To PV"};
  Configurable<double> lambdaSettingcosPA{"lambdaSettingcosPA", 0.98, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> lambdaSettingradius{"lambdaSettingradius", 0.50, "v0radius"};

  // Configurable Antilambda Topological Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> antilambdaSettingdcav0dau{"antilambdaSettingdcav0dau", 0.3, "DCA V0 Daughters"};
  Configurable<float> antilambdaSettingdcapostopv{"antilambdaSettingdcapostopv", 0.09, "DCA Pos To PV"};
  Configurable<float> antilambdaSettingdcanegtopv{"antilambdaSettingdcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<double> antilambdaSettingcosPA{"antilambdaSettingcosPA", 0.98, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> antilambdaSettingradius{"antilambdaSettingradius", 0.50, "v0radius"};

  // Configurables for Specific V0s analysis
  Configurable<bool> kzeroAnalysis{"kzeroAnalysis", true, "Enable Kzerosh Pt Analysis"};
  Configurable<bool> lambdaAnalysis{"lambdaAnalysis", true, "Enable Lambda Pt Analysis"};
  Configurable<bool> antiLambdaAnalysis{"antiLambdaAnalysis", true, "Enable Antilambda Pt Analysis"};

  // Configurable string for Different Pt Bins
  Configurable<std::string> kzeroSettingPtBinsString{"kzeroSettingPtBinsString", {"0.0,0.15,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85,3.0"}, "Kzero Pt Bin Values"};
  Configurable<std::string> lambdaSettingPtBinsString{"lambdaSettingPtBinsString", {"0.0,0.15,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85,3.0"}, "Lambda Pt Bin Values"};
  Configurable<std::string> antilambdaSettingPtBinsString{"antilambdaSettingPtBinsString", {"0.0,0.15,0.3,0.45,0.6,0.75,0.9,1.05,1.2,1.35,1.5,1.65,1.8,1.95,2.1,2.25,2.4,2.55,2.7,2.85,3.0"}, "Antilambda Pt Bin Values"};

  void init(InitContext const&)
  {
    pthistos::kaonPt.resize(nmaxHistograms);       // number of Kaon Pt histograms to expect
    pthistos::lambdaPt.resize(nmaxHistograms);     // number of Lambda histograms to expect
    pthistos::antilambdaPt.resize(nmaxHistograms); // number of Antilambda histograms to expect
    // tokenise strings into individual values
    pthistos::kaonPtBins = o2::utils::Str::tokenize(kzeroSettingPtBinsString, ',');
    pthistos::lambdaPtBins = o2::utils::Str::tokenize(lambdaSettingPtBinsString, ',');
    pthistos::antilambdaPtBins = o2::utils::Str::tokenize(antilambdaSettingPtBinsString, ',');

    // initialize and convert tokenized strings into vector of doubles for AxisSpec
    std::vector<double> kaonptedgevalues(nmaxHistograms + 1);
    std::vector<double> lambdaptedgevalues(nmaxHistograms + 1);
    std::vector<double> antilambdaPtedgevalues(nmaxHistograms + 1);
    for (int i = 0; i < nmaxHistograms + 1; i++) {
      kaonptedgevalues[i] = std::stod(pthistos::kaonPtBins[i]);
      lambdaptedgevalues[i] = std::stod(pthistos::lambdaPtBins[i]);
      antilambdaPtedgevalues[i] = std::stod(pthistos::antilambdaPtBins[i]);
    }

    // Axes
    AxisSpec k0ShortMassAxis = {nBins, 0.45f, 0.55f, "#it{M} #pi^{+}#pi^{-} [GeV/#it{c}^{2}]"};
    AxisSpec lambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M} p^{+}#pi^{-} [GeV/#it{c}^{2}]"};
    AxisSpec antiLambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M} p^{-}#pi^{+} [GeV/#it{c}^{2}]"};
    AxisSpec k0ShortPtAxis = {kaonptedgevalues, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec lambdaPtAxis = {lambdaptedgevalues, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec antilambdaPtAxis = {antilambdaPtedgevalues, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec armenterosQtAxis = {nBinsArmenteros, 0.0f, 0.3f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec armenterosasymAxis = {nBinsArmenteros, -1.f, 1.f, "#ait{p}^{+}_{||}-it{p}^{-}_{||}/it{p}^{+}_{||}+it{p}^{-}_{||}"};
    AxisSpec vertexZAxis = {nBins, -10.0f, 10.0f, "vrtx_{Z} [cm]"};

    std::vector<std::string> kaonhistvalue(nmaxHistograms + 1);
    std::vector<std::string> lambdahistvalue(nmaxHistograms + 1);
    std::vector<std::string> antilambdahistvalue(nmaxHistograms + 1);
    // K0short Histogram Pt Bin Edges
    for (int i = 0; i < nmaxHistograms + 1; i++) {     // Histos won't accept "." character so converting it to "_"
      std::string kaonptbin = pthistos::kaonPtBins[i]; // getting the value of the bin edge
      size_t pos = kaonptbin.find(".");                // finding the "." character
      kaonptbin[pos] = '_';                            // changing the "." character of thestring-value to a "_"
      kaonhistvalue[i] = kaonptbin;                    // filling bin edges list
    }
    // Lambda Histograms Pt Bin Edges (same as K0s above)
    for (int i = 0; i < nmaxHistograms + 1; i++) {
      std::string lambdaptbin = pthistos::lambdaPtBins[i];
      size_t pos = lambdaptbin.find(".");
      lambdaptbin[pos] = '_';
      lambdahistvalue[i] = lambdaptbin;
    }
    // AntiLambda Histograms Pt Bin Edges (same as K0s above)
    for (int i = 0; i < nmaxHistograms + 1; i++) {
      std::string antilambdaPtbin = pthistos::antilambdaPtBins[i];
      size_t pos = antilambdaPtbin.find(".");
      antilambdaPtbin[pos] = '_';
      antilambdahistvalue[i] = antilambdaPtbin;
    }

    rPtAnalysis.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rPtAnalysis.add("hArmenterosPodolanskiPlot", "hArmenterosPodolanskiPlot", {HistType::kTH2F, {{armenterosasymAxis}, {armenterosQtAxis}}});
    rPtAnalysis.add("hV0EtaDaughters", "hV0EtaDaughters", {HistType::kTH1F, {{nBins, -1.2f, 1.2f}}});

    // Generated Pt Spectrums For Feeddown
    rPtAnalysis.add("hXiMinusGeneratedPtSpectrum", "hXiMinusGeneratedPtSpectrum", {HistType::kTH1F, {lambdaPtAxis}});
    rPtAnalysis.add("hXiZeroGeneratedPtSpectrum", "hXiZeroGeneratedPtSpectrum", {HistType::kTH1F, {lambdaPtAxis}});
    rPtAnalysis.add("hOmegaGeneratedPtSpectrum", "hOmegaGeneratedPtSpectrum", {HistType::kTH1F, {lambdaPtAxis}});
    rPtAnalysis.add("hXiPlusGeneratedPtSpectrum", "hXiPlusGeneratedPtSpectrum", {HistType::kTH1F, {antilambdaPtAxis}});
    rPtAnalysis.add("hAntiXiZeroGeneratedPtSpectrum", "hAntiXiZeroGeneratedPtSpectrum", {HistType::kTH1F, {antilambdaPtAxis}});
    rPtAnalysis.add("hAntiOmegaGeneratedPtSpectrum", "hAntiOmegaGeneratedPtSpectrum", {HistType::kTH1F, {antilambdaPtAxis}});

    // Adding Kzerosh Histograms to registry
    if (kzeroAnalysis == true) {
      rPtAnalysis.add("hK0ShGeneratedPtSpectrum", "hK0ShGeneratedPtSpectrum", {HistType::kTH1F, {k0ShortPtAxis}});
      rPtAnalysis.add("hK0ShortReconstructedPtSpectrum", "hK0ShortReconstructedPtSpectrum", {HistType::kTH1F, {k0ShortPtAxis}});
      rPtAnalysis.add("hMassK0ShortAll", "hMassK0ShortAll", {HistType::kTH1F, {k0ShortMassAxis}});
      rPtAnalysis.add("hK0ShortPtSpectrumBeforeCuts", "hK0ShortPtSpectrumBeforeCuts", {HistType::kTH1F, {k0ShortPtAxis}});
      rPtAnalysis.add("hMassK0ShortAllAfterCuts", "hMassK0ShortAllAfterCuts", {HistType::kTH1F, {k0ShortMassAxis}});
      rPtAnalysis.add("hNSigmaPosPiFromK0s", "hNSigmaPosPiFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {k0ShortPtAxis}}});
      rPtAnalysis.add("hNSigmaNegPiFromK0s", "hNSigmaNegPiFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {k0ShortPtAxis}}});
      rPtAnalysis.add("hK0shEtaDaughters", "hK0shEtaDaughters", {HistType::kTH1F, {{nBins, -1.2f, 1.2f}}});
      rPtAnalysis.add("hArmenterosPodolanskiPlotK0Short", "hArmenterosPodolanskiPlotK0Short", {HistType::kTH2F, {{armenterosasymAxis}, {armenterosQtAxis}}});
      for (int i = 0; i < nmaxHistograms; i++) {
        pthistos::kaonPt[i] = rKaonshMassPlotsPerPtBin.add<TH1>(fmt::format("hPt_from_{0}_to_{1}", kaonhistvalue[i], kaonhistvalue[i + 1]).c_str(), fmt::format("hPt_from_{0}_to_{1}", kaonhistvalue[i], kaonhistvalue[i + 1]).c_str(), {HistType::kTH1D, {{k0ShortMassAxis}}});
      }
    }
    // Adding Lambda Histograms
    if (lambdaAnalysis == true) {
      // same method as in Kzerosh above
      rPtAnalysis.add("hLambdaGeneratedPtSpectrum", "hLambdaGeneratedPtSpectrum", {HistType::kTH1F, {lambdaPtAxis}});
      rPtAnalysis.add("hLambdaReconstructedPtSpectrum", "hLambdaReconstructedPtSpectrum", {HistType::kTH1F, {lambdaPtAxis}});
      rPtAnalysis.add("hMassLambdaAll", "hMassLambdaAll", {HistType::kTH1F, {lambdaMassAxis}});
      rPtAnalysis.add("hLambdaPtSpectrumBeforeCuts", "hLambdaPtSpectrumBeforeCuts", {HistType::kTH1F, {lambdaPtAxis}});
      rPtAnalysis.add("hMassLambdaAllAfterCuts", "hMassLambdaAllAfterCuts", {HistType::kTH1F, {lambdaMassAxis}});
      rPtAnalysis.add("hNSigmaPosProtonFromLambda", "hNSigmaPosProtonFromLambda", {HistType::kTH2F, {{100, -5.f, 5.f}, {lambdaPtAxis}}});
      rPtAnalysis.add("hNSigmaNegPionFromLambda", "hNSigmaNegPionFromLambda", {HistType::kTH2F, {{100, -5.f, 5.f}, {lambdaPtAxis}}});
      rPtAnalysis.add("hLambdaEtaDaughters", "hLambdaEtaDaughters", {HistType::kTH1F, {{nBins, -1.2f, 1.2f}}});
      rPtAnalysis.add("hArmenterosPodolanskiPlotLambda", "hArmenterosPodolanskiPlotLambda", {HistType::kTH2F, {{armenterosasymAxis}, {armenterosQtAxis}}});
      for (int i = 0; i < nmaxHistograms; i++) {
        pthistos::lambdaPt[i] = rLambdaMassPlotsPerPtBin.add<TH1>(fmt::format("hPt_from_{0}_to_{1}", lambdahistvalue[i], lambdahistvalue[i + 1]).c_str(), fmt::format("hPt_from_{0}_to_{1}", lambdahistvalue[i], lambdahistvalue[i + 1]).c_str(), {HistType::kTH1D, {{lambdaMassAxis}}});
      }
      // lambdafeeddown matrices
      rFeeddownMatrices.add("hLambdaXiMinusFeeddownMatrix", "hLambdaXiMinusFeeddownMatrix", {HistType::kTH2F, {{lambdaPtAxis}, {lambdaPtAxis}}});
      rFeeddownMatrices.add("hLambdaXiZeroFeeddownMatrix", "hLambdaXiZeroFeeddownMatrix", {HistType::kTH2F, {{lambdaPtAxis}, {lambdaPtAxis}}});
      rFeeddownMatrices.add("hLambdaOmegaFeeddownMatrix", "hLambdaOmegaFeeddownMatrix", {HistType::kTH2F, {{lambdaPtAxis}, {lambdaPtAxis}}});
    }
    // Adding Antilambda Histograms
    if (antiLambdaAnalysis == true) {
      // same method as in Lambda and Kzerosh above
      rPtAnalysis.add("hAntilambdaGeneratedPtSpectrum", "hAntilambdaGeneratedPtSpectrum", {HistType::kTH1F, {{antilambdaPtAxis}}});
      rPtAnalysis.add("hAntilambdaReconstructedPtSpectrum", "hAntilambdaReconstructedPtSpectrum", {HistType::kTH1F, {antilambdaPtAxis}});
      rPtAnalysis.add("hMassAntilambdaAll", "hMassAntilambdaAll", {HistType::kTH1F, {antiLambdaMassAxis}});
      rPtAnalysis.add("hantilambdaPtSpectrumBeforeCuts", "hantilambdaPtSpectrumBeforeCuts", {HistType::kTH1F, {antilambdaPtAxis}});
      rPtAnalysis.add("hMassAntilambdaAllAfterCuts", "hMassAntilambdaAllAfterCuts", {HistType::kTH1F, {antiLambdaMassAxis}});
      rPtAnalysis.add("hNSigmaNegProtonFromAntilambda", "hNSigmaNegProtonFromAntilambda", {HistType::kTH2F, {{100, -5.f, 5.f}, {antilambdaPtAxis}}});
      rPtAnalysis.add("hNSigmaPosPionFromAntilambda", "hNSigmaPosPionFromAntilambda", {HistType::kTH2F, {{100, -5.f, 5.f}, {antilambdaPtAxis}}});
      rPtAnalysis.add("hAntiLambdaEtaDaughters", "hAntiLambdaEtaDaughters", {HistType::kTH1F, {{nBins, -1.2f, 1.2f}}});
      rPtAnalysis.add("hArmenterosPodolanskiPlotAntiLambda", "hArmenterosPodolanskiPlotAntiLambda", {HistType::kTH2F, {{armenterosasymAxis}, {armenterosQtAxis}}});
      for (int i = 0; i < nmaxHistograms; i++) {
        pthistos::antilambdaPt[i] = rAntilambdaMassPlotsPerPtBin.add<TH1>(fmt::format("hPt_from_{0}_to_{1}", antilambdahistvalue[i], antilambdahistvalue[i + 1]).c_str(), fmt::format("hPt_from_{0}_to_{1}", antilambdahistvalue[i], antilambdahistvalue[i + 1]).c_str(), {HistType::kTH1D, {{antiLambdaMassAxis}}});
      }
      // antilambdafeeddown matrices
      rFeeddownMatrices.add("hAntiLambdaXiPlusFeeddownMatrix", "hAntiLambdaXiPlusFeeddownMatrix", {HistType::kTH2F, {{antilambdaPtAxis}, {antilambdaPtAxis}}});
      rFeeddownMatrices.add("hAntiLambdaAntiXiZeroFeeddownMatrix", "hAntiLambdaAntiXiZeroFeeddownMatrix", {HistType::kTH2F, {{antilambdaPtAxis}, {antilambdaPtAxis}}});
      rFeeddownMatrices.add("hAntiLambdaAntiOmegaFeeddownMatrix", "hAntiLambdaAntiOmegaPlusFeeddownMatrix", {HistType::kTH2F, {{antilambdaPtAxis}, {antilambdaPtAxis}}});
    }
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < cutZVertex);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutZVertex);

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;

  // This is the Process for the MC Generated Data
  void genMCProcess(soa::Filtered<aod::McCollisions>::iterator const&,
                    const soa::SmallGroups<soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>>&,
                    aod::McParticles const& mcParticles)
  {
    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) < rapiditycutGen) {
        if (mcParticle.isPhysicalPrimary()) {
          if (mcParticle.pdgCode() == 310) // kzero matched
          {
            rPtAnalysis.fill(HIST("hK0ShGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == 3122) // lambda matched
          {
            rPtAnalysis.fill(HIST("hLambdaGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == -3122) // antilambda matched
          {
            rPtAnalysis.fill(HIST("hAntilambdaGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == 3312) // Xi Minus matched
          {
            rPtAnalysis.fill(HIST("hXiMinusGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == 3322) // Xi Zero matched
          {
            rPtAnalysis.fill(HIST("hXiZeroGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == 3334) // Omega matched
          {
            rPtAnalysis.fill(HIST("hOmegaGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == -3312) // Xi Plus matched
          {
            rPtAnalysis.fill(HIST("hXiPlusGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == -3322) // Anti-Xi Zero matched
          {
            rPtAnalysis.fill(HIST("hAntiXiZeroGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == -3334) // Anti-Omega matched
          {
            rPtAnalysis.fill(HIST("hAntiOmegaGeneratedPtSpectrum"), mcParticle.pt());
          }
        }
      }
    }
  }
  // This is the Process for the MC reconstructed Data
  void recMCProcess(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                    soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s,
                    DaughterTracks const&, // no need to define a variable for tracks, if we don't access them directly
                    aod::McParticles const&)
  {
    // PDG mass values for Competitive V0 Cut old: const auto& mK0shPDG = 0.497611;
    double mK0shPDG = o2::constants::physics::MassK0Short;
    double mLambdaPDG = o2::constants::physics::MassLambda0;

    // tokenise strings into individual values
    pthistos::kaonPtBins = o2::utils::Str::tokenize(kzeroSettingPtBinsString, ',');
    pthistos::lambdaPtBins = o2::utils::Str::tokenize(lambdaSettingPtBinsString, ',');
    pthistos::antilambdaPtBins = o2::utils::Str::tokenize(antilambdaSettingPtBinsString, ',');

    // initialize and convert tokenized strings into vector of doubles for Pt Bin Edges
    std::vector<double> kaonptedgevalues(nmaxHistograms + 1);
    std::vector<double> lambdaptedgevalues(nmaxHistograms + 1);
    std::vector<double> antilambdaPtedgevalues(nmaxHistograms + 1);

    for (int i = 0; i < nmaxHistograms + 1; i++) {
      kaonptedgevalues[i] = std::stod(pthistos::kaonPtBins[i]);
      lambdaptedgevalues[i] = std::stod(pthistos::lambdaPtBins[i]);
      antilambdaPtedgevalues[i] = std::stod(pthistos::antilambdaPtBins[i]);
    }

    for (const auto& v0 : V0s) {
      rPtAnalysis.fill(HIST("hVertexZRec"), collision.posZ());
      // Checking that the V0 is a true K0s/Lambdas/Antilambdas and then filling the parameter histograms and the invariant mass plots for different cuts (which are taken from namespace)
      if (v0.has_mcParticle()) {
        auto v0mcParticle = v0.mcParticle();
        if (std::abs(v0.posTrack_as<DaughterTracks>().eta()) < etadau && std::abs(v0.negTrack_as<DaughterTracks>().eta()) < etadau) { // daughters pseudorapidity cut
          rPtAnalysis.fill(HIST("hV0EtaDaughters"), v0.negTrack_as<DaughterTracks>().eta());
          rPtAnalysis.fill(HIST("hV0EtaDaughters"), v0.posTrack_as<DaughterTracks>().eta());
          if (kzeroAnalysis == true) {
            if (v0mcParticle.pdgCode() == 310) { // kzero matched
              rPtAnalysis.fill(HIST("hMassK0ShortAll"), v0.mK0Short());
              rPtAnalysis.fill(HIST("hK0ShortPtSpectrumBeforeCuts"), v0.pt());
              if (std::abs(v0.mLambda() - mLambdaPDG) > compv0masscut && std::abs(v0.mAntiLambda() - mLambdaPDG) > compv0masscut) { // Kzero competitive v0 mass cut (cut out Lambdas and Anti-Lambdas)
                // Implementing best kzero topological cuts
                if (v0.v0cosPA() > kaonshSettingcosPA && v0.dcaV0daughters() < kaonshSettingdcav0dau && v0.v0radius() > kaonshSettingradius && std::abs(v0.dcapostopv()) > kaonshSettingdcapostopv && std::abs(v0.dcanegtopv()) > kaonshSettingdcanegtopv) {
                  rPtAnalysis.fill(HIST("hMassK0ShortAllAfterCuts"), v0.mK0Short());
                  rPtAnalysis.fill(HIST("hK0ShortReconstructedPtSpectrum"), v0.pt());
                  rPtAnalysis.fill(HIST("hK0shEtaDaughters"), v0.negTrack_as<DaughterTracks>().eta());
                  rPtAnalysis.fill(HIST("hK0shEtaDaughters"), v0.posTrack_as<DaughterTracks>().eta());
                  if (v0mcParticle.isPhysicalPrimary()) {
                    for (int i = 0; i < nmaxHistograms; i++) {
                      if (kaonptedgevalues[i] <= v0.pt() && v0.pt() < kaonptedgevalues[i + 1]) { // finding v0s with pt within the range of our bin edges
                        pthistos::kaonPt[i]->Fill(v0.mK0Short());                                // filling the k0s namespace histograms
                      }
                    }
                  }
                }
              }
            }
          }
          // lambda analysis
          if (lambdaAnalysis == true) {
            if (v0mcParticle.pdgCode() == 3122) { // lambda matched
              rPtAnalysis.fill(HIST("hMassLambdaAll"), v0.mLambda());
              rPtAnalysis.fill(HIST("hLambdaPtSpectrumBeforeCuts"), v0.pt());
              if (std::abs(v0.mK0Short() - mK0shPDG) > compv0masscut) { // antilambda competitive v0 mass cut (cut out Kaons)
                // Implementing best lambda cuts
                if (v0.v0cosPA() > lambdaSettingcosPA && v0.dcaV0daughters() < lambdaSettingdcav0dau && v0.v0radius() > lambdaSettingradius && std::abs(v0.dcapostopv()) > lambdaSettingdcapostopv && std::abs(v0.dcanegtopv()) > lambdaSettingdcanegtopv) {
                  rPtAnalysis.fill(HIST("hMassLambdaAllAfterCuts"), v0.mLambda());
                  rPtAnalysis.fill(HIST("hLambdaReconstructedPtSpectrum"), v0.pt());
                  rPtAnalysis.fill(HIST("hLambdaEtaDaughters"), v0.negTrack_as<DaughterTracks>().eta());
                  rPtAnalysis.fill(HIST("hLambdaEtaDaughters"), v0.posTrack_as<DaughterTracks>().eta());
                  if (v0mcParticle.isPhysicalPrimary()) {
                    for (int i = 0; i < nmaxHistograms; i++) {
                      if (lambdaptedgevalues[i] <= v0.pt() && v0.pt() < lambdaptedgevalues[i + 1]) {
                        pthistos::lambdaPt[i]->Fill(v0.mLambda());
                      }
                    }
                  }
                  if (!v0mcParticle.isPhysicalPrimary()) {
                    auto v0mothers = v0mcParticle.mothers_as<aod::McParticles>(); // Get mothers
                    if (!v0mothers.empty()) {
                      auto& v0mcParticleMother = v0mothers.front(); // First mother
                      if (v0mcParticleMother.pdgCode() == 3312)     // Xi Minus Mother Matched
                      {
                        rFeeddownMatrices.fill(HIST("hLambdaXiMinusFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                      }
                      if (v0mcParticleMother.pdgCode() == 3322) // Xi Zero Mother Matched
                      {
                        rFeeddownMatrices.fill(HIST("hLambdaXiZeroFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                      }
                      if (v0mcParticleMother.pdgCode() == 3334) // Omega Mother Matched
                      {
                        rFeeddownMatrices.fill(HIST("hLambdaOmegaFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                      }
                    }
                  }
                }
              }
            }
          }
          // antilambda analysis
          if (antiLambdaAnalysis == true) {
            if (v0mcParticle.pdgCode() == -3122) { // antilambda matched
              rPtAnalysis.fill(HIST("hMassAntilambdaAll"), v0.mAntiLambda());
              rPtAnalysis.fill(HIST("hantilambdaPtSpectrumBeforeCuts"), v0.pt());
              if (std::abs(v0.mK0Short() - mK0shPDG) > compv0masscut) { // antilambda competitive v0 mass cut (cut out Kaons)
                // Implementing best antilambda cuts
                if (v0.v0cosPA() > antilambdaSettingcosPA && v0.dcaV0daughters() < antilambdaSettingdcav0dau && v0.v0radius() > antilambdaSettingradius && std::abs(v0.dcapostopv()) > antilambdaSettingdcapostopv && std::abs(v0.dcanegtopv()) > antilambdaSettingdcanegtopv) {
                  rPtAnalysis.fill(HIST("hMassAntilambdaAllAfterCuts"), v0.mAntiLambda());
                  rPtAnalysis.fill(HIST("hAntilambdaReconstructedPtSpectrum"), v0.pt());
                  rPtAnalysis.fill(HIST("hAntiLambdaEtaDaughters"), v0.negTrack_as<DaughterTracks>().eta());
                  rPtAnalysis.fill(HIST("hAntiLambdaEtaDaughters"), v0.posTrack_as<DaughterTracks>().eta());
                  if (v0mcParticle.isPhysicalPrimary()) {
                    for (int i = 0; i < nmaxHistograms; i++) {
                      if (antilambdaPtedgevalues[i] <= v0.pt() && v0.pt() < antilambdaPtedgevalues[i + 1]) {
                        pthistos::antilambdaPt[i]->Fill(v0.mAntiLambda());
                      }
                    }
                  }
                  if (!v0mcParticle.isPhysicalPrimary()) {
                    auto v0mothers = v0mcParticle.mothers_as<aod::McParticles>(); // Get mothers
                    if (!v0mothers.empty()) {
                      auto& v0mcParticleMother = v0mothers.front(); // First mother
                      if (v0mcParticleMother.pdgCode() == -3312)    // Xi Plus Mother Matched
                      {
                        rFeeddownMatrices.fill(HIST("hAntiLambdaXiPlusFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                      }
                      if (v0mcParticleMother.pdgCode() == -3322) // Anti-Xi Zero Mother Matched
                      {
                        rFeeddownMatrices.fill(HIST("hAntiLambdaAntiXiZeroFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                      }
                      if (v0mcParticleMother.pdgCode() == -3334) // Anti-Omega (minus) Mother Matched
                      {
                        rFeeddownMatrices.fill(HIST("hAntiLambdaAntiOmegaFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  // This is the process for Real Data
  void dataProcess(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                   aod::V0Datas const& V0s,
                   DaughterTracks const&)
  {
    double mK0shPDG = o2::constants::physics::MassK0Short;
    double mLambdaPDG = o2::constants::physics::MassLambda0;

    // tokenise strings into individual values
    pthistos::kaonPtBins = o2::utils::Str::tokenize(kzeroSettingPtBinsString, ',');
    pthistos::lambdaPtBins = o2::utils::Str::tokenize(lambdaSettingPtBinsString, ',');
    pthistos::antilambdaPtBins = o2::utils::Str::tokenize(antilambdaSettingPtBinsString, ',');

    // initialize and convert tokenized strings into vector of doubles for pt bin edges
    std::vector<double> kaonptedgevalues(nmaxHistograms + 1);
    std::vector<double> lambdaptedgevalues(nmaxHistograms + 1);
    std::vector<double> antilambdaPtedgevalues(nmaxHistograms + 1);
    for (int i = 0; i < nmaxHistograms + 1; i++) {
      kaonptedgevalues[i] = std::stod(pthistos::kaonPtBins[i]);
      lambdaptedgevalues[i] = std::stod(pthistos::lambdaPtBins[i]);
      antilambdaPtedgevalues[i] = std::stod(pthistos::antilambdaPtBins[i]);
    }

    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<DaughterTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<DaughterTracks>();
      // Armenteros-Podolandski Plot Values
      double pv0 = std::sqrt((v0.px() * v0.px()) + (v0.py() * v0.py()) + (v0.pz() * v0.pz()));
      double pposdauparallelv0 = ((v0.posTrack_as<DaughterTracks>().px() * v0.px()) + (v0.posTrack_as<DaughterTracks>().py() * v0.py()) + (v0.posTrack_as<DaughterTracks>().pz() * v0.pz())) / pv0;
      double qValue = std::sqrt(((v0.posTrack_as<DaughterTracks>().px() * v0.posTrack_as<DaughterTracks>().px()) + (v0.posTrack_as<DaughterTracks>().py() * v0.posTrack_as<DaughterTracks>().py()) + (v0.posTrack_as<DaughterTracks>().pz() * v0.posTrack_as<DaughterTracks>().pz())) - (pposdauparallelv0 * pposdauparallelv0));
      double plpos = (v0.posTrack_as<DaughterTracks>().px() * v0.px() / pv0) + (v0.posTrack_as<DaughterTracks>().py() * v0.py() / pv0) + (v0.posTrack_as<DaughterTracks>().pz() * v0.pz() / pv0);
      double plneg = (v0.negTrack_as<DaughterTracks>().px() * v0.px() / pv0) + (v0.negTrack_as<DaughterTracks>().py() * v0.py() / pv0) + (v0.negTrack_as<DaughterTracks>().pz() * v0.pz() / pv0);
      double aValue = (plpos - plneg) / (plpos + plneg);
      rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlot"), aValue, qValue);
      rPtAnalysis.fill(HIST("hVertexZRec"), collision.posZ());
      if (std::abs(v0.posTrack_as<DaughterTracks>().eta()) < etadau && std::abs(v0.negTrack_as<DaughterTracks>().eta()) < etadau) { // daughters pseudorapidity cut
        rPtAnalysis.fill(HIST("hV0EtaDaughters"), v0.negTrack_as<DaughterTracks>().eta());
        rPtAnalysis.fill(HIST("hV0EtaDaughters"), v0.posTrack_as<DaughterTracks>().eta());
        // kzero analysis
        if (kzeroAnalysis == true) {
          // Filling the five Kzero invariant mass plots for different cuts (which are taken from namespace), for full explanation see the first kzero cut filling in the MC process
          rPtAnalysis.fill(HIST("hMassK0ShortAll"), v0.mK0Short());
          if (std::abs(v0.mLambda() - mLambdaPDG) > compv0masscut && std::abs(v0.mAntiLambda() - mLambdaPDG) > compv0masscut) { // antilambda competitive v0 mass cut (cut out Lambdas and Anti-Lambdas)
            // Implementing best kzero cuts
            if (std::abs(posDaughterTrack.tpcNSigmaPi()) < nSigmaTPCPion && std::abs(negDaughterTrack.tpcNSigmaPi()) < nSigmaTPCPion) { // TPC PID on daughter pions
              rPtAnalysis.fill(HIST("hNSigmaPosPiFromK0s"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tpcInnerParam());
              rPtAnalysis.fill(HIST("hNSigmaNegPiFromK0s"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
              if (v0.v0cosPA() > kaonshSettingcosPA && v0.dcaV0daughters() < kaonshSettingdcav0dau && v0.v0radius() > kaonshSettingradius && std::abs(v0.dcapostopv()) > kaonshSettingdcapostopv && std::abs(v0.dcanegtopv()) > kaonshSettingdcanegtopv) {
                rPtAnalysis.fill(HIST("hMassK0ShortAllAfterCuts"), v0.mK0Short());
                rPtAnalysis.fill(HIST("hK0shEtaDaughters"), v0.negTrack_as<DaughterTracks>().eta());
                rPtAnalysis.fill(HIST("hK0shEtaDaughters"), v0.posTrack_as<DaughterTracks>().eta());
                rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlotK0Short"), aValue, qValue);
                for (int i = 0; i < nmaxHistograms; i++) {
                  if (kaonptedgevalues[i] <= v0.pt() && v0.pt() < kaonptedgevalues[i + 1]) {
                    pthistos::kaonPt[i]->Fill(v0.mK0Short());
                  }
                }
              }
            }
          }
        }
        // lambda analysis
        if (lambdaAnalysis == true) {
          // Filling the five lambda invariant mass plots for different cuts (which are taken from namespace), for full explanation see the first kzero cut filling in the MC process
          rPtAnalysis.fill(HIST("hMassLambdaAll"), v0.mLambda());
          if (std::abs(v0.mK0Short() - mK0shPDG) > compv0masscut) {                                                                       // lambda competitive v0 mass cut (cut out Kaons)
            if (std::abs(posDaughterTrack.tpcNSigmaPr()) < nSigmaTPCProton && std::abs(negDaughterTrack.tpcNSigmaPi()) < nSigmaTPCPion) { // TPC PID on daughter pion and proton for Lambda
              rPtAnalysis.fill(HIST("hNSigmaPosProtonFromLambda"), posDaughterTrack.tpcNSigmaPr(), posDaughterTrack.tpcInnerParam());
              rPtAnalysis.fill(HIST("hNSigmaNegPionFromLambda"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
              // Implementing best lambda cuts
              if (v0.v0cosPA() > lambdaSettingcosPA && v0.dcaV0daughters() < lambdaSettingdcav0dau && v0.v0radius() > lambdaSettingradius && std::abs(v0.dcapostopv()) > lambdaSettingdcapostopv && std::abs(v0.dcanegtopv()) > lambdaSettingdcanegtopv) {
                rPtAnalysis.fill(HIST("hMassLambdaAllAfterCuts"), v0.mLambda());
                rPtAnalysis.fill(HIST("hLambdaEtaDaughters"), v0.negTrack_as<DaughterTracks>().eta());
                rPtAnalysis.fill(HIST("hLambdaEtaDaughters"), v0.posTrack_as<DaughterTracks>().eta());
                rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlotLambda"), aValue, qValue);
                for (int i = 0; i < nmaxHistograms; i++) {
                  if (lambdaptedgevalues[i] <= v0.pt() && v0.pt() < lambdaptedgevalues[i + 1]) {
                    pthistos::lambdaPt[i]->Fill(v0.mLambda());
                  }
                }
              }
            }
          }
        }
        // anti-lambda analysis
        if (antiLambdaAnalysis == true) {
          // Filling the five Antilambda invariant mass plots for different cuts (which are taken from namespace), for full explanation see the first kzero cut filling in the MC process
          rPtAnalysis.fill(HIST("hMassAntilambdaAll"), v0.mAntiLambda());
          if (std::abs(v0.mK0Short() - mK0shPDG) > compv0masscut) {                                                                       // antilambda competitive v0 mass cut (cut out Kaons)
            if (std::abs(negDaughterTrack.tpcNSigmaPr()) < nSigmaTPCProton && std::abs(posDaughterTrack.tpcNSigmaPi()) < nSigmaTPCPion) { // TPC PID on daughter pion and proton for AntiLambda
              rPtAnalysis.fill(HIST("hNSigmaPosPionFromAntilambda"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tpcInnerParam());
              rPtAnalysis.fill(HIST("hNSigmaNegProtonFromAntilambda"), negDaughterTrack.tpcNSigmaPr(), negDaughterTrack.tpcInnerParam());
              // implementing best antilambda cuts
              if (v0.v0cosPA() > antilambdaSettingcosPA && v0.dcaV0daughters() < antilambdaSettingdcav0dau && v0.v0radius() > antilambdaSettingradius && std::abs(v0.dcapostopv()) > antilambdaSettingdcapostopv && std::abs(v0.dcanegtopv()) > antilambdaSettingdcanegtopv) {
                rPtAnalysis.fill(HIST("hMassAntilambdaAllAfterCuts"), v0.mAntiLambda());
                rPtAnalysis.fill(HIST("hAntiLambdaEtaDaughters"), v0.negTrack_as<DaughterTracks>().eta());
                rPtAnalysis.fill(HIST("hAntiLambdaEtaDaughters"), v0.posTrack_as<DaughterTracks>().eta());
                rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlotAntiLambda"), aValue, qValue);
                for (int i = 0; i < nmaxHistograms; i++) {
                  if (lambdaptedgevalues[i] <= v0.pt() && v0.pt() < lambdaptedgevalues[i + 1]) {
                    pthistos::antilambdaPt[i]->Fill(v0.mAntiLambda());
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(V0PtInvMassPlots, genMCProcess, "Process Run 3 MC Generated", false);
  PROCESS_SWITCH(V0PtInvMassPlots, recMCProcess, "Process Run 3 MC Reconstructed", false);
  PROCESS_SWITCH(V0PtInvMassPlots, dataProcess, "Process Run 3 Data,", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<V0PtInvMassPlots>(cfgc)};
}
