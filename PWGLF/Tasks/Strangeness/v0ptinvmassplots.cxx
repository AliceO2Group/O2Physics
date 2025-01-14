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
This task creates 20 histograms that are filled with the V0 invariant mass under the K0, Lambda and Antilambda mass assumption
for different pt ranges (constituting bins), so 3x20=60 plots.The values are inserted as configurable strings for convinience.
Plots of the invariant masses at different stages of the analysis (ex. before and after the V0 cuts are enforced) and some pt distributions.
This analysis includes two processes, one for Real Data and one for MC Data switchable at the end of the code, only run one at a time*/

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
constexpr uint32_t NSIZE = 30;
std::shared_ptr<TH1> kaonPt[NSIZE];
static std::vector<std::string> kaonptbins;
std::shared_ptr<TH1> lambdaPt[NSIZE];
static std::vector<std::string> lambdaPtBins;
std::shared_ptr<TH1> antiLambdaPt[NSIZE];
static std::vector<std::string> antiLambdaPtBins;
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

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  Configurable<int> xaxisGenBins{"xaxisGenBins", pthistos::NSIZE, "Number of bins for Generated Pt Spectrum"}; // is nSIZE ok
  Configurable<float> xaxisMinGenBin{"xaxisMinGenBin", 0.0, "Minimum bin value of the Generated Pt Spectrum Plot"};
  Configurable<float> xaxisMaxGenBin{"xaxisMaxGenBin", 3.0, "Maximum bin value of the Generated Pt Spectrum Plot"};

  // Configurable Kaonsh Topological Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> kaonshSettingdcav0dau{"kaonshSettingdcav0dau", 100.0, "DCA V0 Daughters"};
  Configurable<float> kaonshSettingdcapostopv{"kaonshSettingdcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<float> kaonshSettingdcanegtopv{"kaonshSettingdcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<double> kaonshSettingcosPA{"kaonshSettingcosPA", 0.50, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> kaonshSettingradius{"kaonshSettingradius", 0.50, "v0radius"};

  // Configurable Lambda Topological Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> lambdaSettingdcav0dau{"lambdaSettingdcav0dau", 100.0, "DCA V0 Daughters"};
  Configurable<float> lambdaSettingdcapostopv{"lambdaSettingdcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<float> lambdaSettingdcanegtopv{"lambdaSettingdcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<double> lambdaSettingcosPA{"lambdaSettingcosPA", 0.50, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> lambdaSettingradius{"lambdaSettingradius", 0.50, "v0radius"};

  // Configurable Antilambda Topological Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> antilambdaSettingdcav0dau{"antilambdaSettingdcav0dau", 100.0, "DCA V0 Daughters"};
  Configurable<float> antilambdaSettingdcapostopv{"antilambdaSettingdcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<float> antilambdaSettingdcanegtopv{"antilambdaSettingdcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<double> antilambdaSettingcosPA{"antilambdaSettingcosPA", 0.50, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> antilambdaSettingradius{"antilambdaSettingradius", 0.50, "v0radius"};

  // Configurables for Specific V0s analysis
  Configurable<bool> kzeroAnalysis{"kzeroAnalysis", true, "Enable Kzerosh Pt Analysis"};
  Configurable<bool> lambdaAnalysis{"lambdaAnalysis", true, "Enable Lambda Pt Analysis"};
  Configurable<bool> antiLambdaAnalysis{"antiLambdaAnalysis", true, "Enable Antilambda Pt Analysis"};

  // Configurable string for Different Pt Bins
  Configurable<std::string> kzeroshSettingPtString{"kzeroSettingptbins", {"0_0,0_15,0_3,0_45,0_6,0_75,0_9,1_05,1_2,1_35,1_5,1_65,1_8,1_95,2_1,2_25,2_4,2_55,2_7,2_85,3_0"}, "Kzero Pt Bin Values"};
  Configurable<std::string> lambdaSettingPtString{"lambdaSettingptbins", {"0_0,0_15,0_3,0_45,0_6,0_75,0_9,1_05,1_2,1_35,1_5,1_65,1_8,1_95,2_1,2_25,2_4,2_55,2_7,2_85,3_0"}, "Lambda Pt Bin Values"};
  Configurable<std::string> antilambdaSettingPtString{"antilambdaSettingptbins", {"0_0,0_15,0_3,0_45,0_6,0_75,0_9,1_05,1_2,1_35,1_5,1_65,1_8,1_95,2_1,2_25,2_4,2_55,2_7,2_85,3_0"}, "Antilambda Pt Bin Values"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec k0ShortMassAxis = {nBins, 0.45f, 0.55f, "#it{M} #pi^{+}#pi^{-} [GeV/#it{c}^{2}]"};
    AxisSpec lambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M} p^{+}#pi^{-} [GeV/#it{c}^{2}]"};
    AxisSpec antiLambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M} p^{-}#pi^{+} [GeV/#it{c}^{2}]"};
    AxisSpec ptAxis = {nBins, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec genPtAxis = {xaxisGenBins, xaxisMinGenBin, xaxisMaxGenBin, "#it{p}_{T} (GeV/#it{c})"};

    rPtAnalysis.add("hV0PtAll", "hV0PtAll", {HistType::kTH1F, {{nBins, 0.0f, 10.0f}}});

    // Adding Kzerosh Histograms to registry
    if (kzeroAnalysis == true) {
      pthistos::kaonptbins = o2::utils::Str::tokenize(kzeroshSettingPtString, ',');
      rPtAnalysis.add("hK0ShortReconstructedPtSpectrum", "hK0ShortReconstructedPtSpectrum", {HistType::kTH1F, {ptAxis}});
      rPtAnalysis.add("hMassK0ShortAll", "hMassK0ShortAll", {HistType::kTH1F, {k0ShortMassAxis}});
      rPtAnalysis.add("hK0ShortPtSpectrumBeforeCuts", "hK0ShortPtSpectrumBeforeCuts", {HistType::kTH1F, {ptAxis}});
      rPtAnalysis.add("hMassK0ShortAllAfterCuts", "hMassK0ShortAllAfterCuts", {HistType::kTH1F, {k0ShortMassAxis}});
      rPtAnalysis.add("hK0ShGeneratedPtSpectrum", "hK0ShGeneratedPtSpectrum", {HistType::kTH1F, {genPtAxis}});
      rPtAnalysis.add("hLambdaGeneratedPtSpectrum", "hLambdaGeneratedPtSpectrum", {HistType::kTH1F, {genPtAxis}});
      rPtAnalysis.add("hAntilambdaGeneratedPtSpectrum", "hAntilambdaGeneratedPtSpectrum", {HistType::kTH1F, {genPtAxis}});
      for (uint32_t i = 0; i < pthistos::kaonptbins.size() - 1; i++) {
        pthistos::kaonPt[i] = rKaonshMassPlotsPerPtBin.add<TH1>(fmt::format("hPt_from_{0}_to_{1}", pthistos::kaonptbins[i], pthistos::kaonptbins[i + 1]).data(), fmt::format("hPt from {0} to {1}", pthistos::kaonptbins[i], pthistos::kaonptbins[i + 1]).data(), {HistType::kTH1D, {{k0ShortMassAxis}}});
      }
    }
    // Adding Lambda Histograms
    if (lambdaAnalysis == true) {
      // same method as in Kzerosh above
      std::string lambdaSettingptbins = lambdaSettingPtString;
      pthistos::lambdaPtBins = o2::utils::Str::tokenize(lambdaSettingptbins, ',');
      rPtAnalysis.add("hLambdaReconstructedPtSpectrum", "hLambdaReconstructedPtSpectrum", {HistType::kTH1F, {ptAxis}});
      rPtAnalysis.add("hMassLambdaAll", "hMassLambdaAll", {HistType::kTH1F, {lambdaMassAxis}});
      rPtAnalysis.add("hLambdaPtSpectrumBeforeCuts", "hLambdaPtSpectrumBeforeCuts", {HistType::kTH1F, {ptAxis}});
      rPtAnalysis.add("hMassLambdaAllAfterCuts", "hMassLambdaAllAfterCuts", {HistType::kTH1F, {lambdaMassAxis}});
      for (u_int32_t i = 0; i < pthistos::lambdaPtBins.size() - 1; i++) {
        pthistos::lambdaPt[i] = rLambdaMassPlotsPerPtBin.add<TH1>(fmt::format("hPt_from_{0}_to_{1}", pthistos::lambdaPtBins[i], pthistos::lambdaPtBins[i + 1]).data(), fmt::format("hPt from {0} to {1}", pthistos::lambdaPtBins[i], pthistos::lambdaPtBins[i + 1]).data(), {HistType::kTH1D, {{lambdaMassAxis}}});
      }
    }
    // Adding Antilambda Histograms
    if (antiLambdaAnalysis == true) {
      // same method as in Lambda and Kzerosh above
      std::string antilambdaSettingptbins = antilambdaSettingPtString;
      pthistos::antiLambdaPtBins = o2::utils::Str::tokenize(antilambdaSettingptbins, ',');
      rPtAnalysis.add("hAntilambdaReconstructedPtSpectrum", "hAntilambdaReconstructedPtSpectrum", {HistType::kTH1F, {ptAxis}});
      rPtAnalysis.add("hMassAntilambdaAll", "hMassAntilambdaAll", {HistType::kTH1F, {antiLambdaMassAxis}});
      rPtAnalysis.add("hAntilambdaPtSpectrumBeforeCuts", "hAntilambdaPtSpectrumBeforeCuts", {HistType::kTH1F, {ptAxis}});
      rPtAnalysis.add("hMassAntilambdaAllAfterCuts", "hMassAntilambdaAllAfterCuts", {HistType::kTH1F, {antiLambdaMassAxis}});
      for (u_int32_t i = 0; i < pthistos::antiLambdaPtBins.size() - 1; i++) {
        pthistos::antiLambdaPt[i] = rAntilambdaMassPlotsPerPtBin.add<TH1>(fmt::format("hPt_from_{0}_to_{1}", pthistos::antiLambdaPtBins[i], pthistos::antiLambdaPtBins[i + 1]).data(), fmt::format("hPt from {0} to {1}", pthistos::antiLambdaPtBins[i], pthistos::antiLambdaPtBins[i + 1]).data(), {HistType::kTH1D, {{antiLambdaMassAxis}}});
      }
    }
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilterMC = (nabs(o2::aod::mccollision::posZ) < 10.0f);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < 10.0f);

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels>;

  // This is the Process for the MC Generated Data
  void genMCProcess(soa::Filtered<aod::McCollisions>::iterator const&,
                    const soa::SmallGroups<soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>>&,
                    aod::McParticles const& mcParticles)
  {
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary()) {
        if (std::abs(mcParticle.y()) < 0.5f) {
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
        }
      }
    }
  }
  // This is the Process for the MC reconstructed Data
  void recMCProcess(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const&,
                    soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s,
                    DaughterTracks const&, // no need to define a variable for tracks, if we don't access them directly
                    aod::McParticles const& /*mcParticles*/)
  {
    for (const auto& v0 : V0s) {
      rPtAnalysis.fill(HIST("hV0PtAll"), v0.pt());
      // Checking that the V0 is a true K0s/Lambdas/Antilambdas and then filling the parameter histograms and the invariant mass plots for different cuts (which are taken from namespace)
      if (v0.has_mcParticle()) {
        auto v0mcParticle = v0.mcParticle();
        if (v0mcParticle.isPhysicalPrimary()) {
          if (kzeroAnalysis == true) {
            if (v0mcParticle.pdgCode() == 310) { // kzero matched
              rPtAnalysis.fill(HIST("hMassK0ShortAll"), v0.mK0Short());
              rPtAnalysis.fill(HIST("hK0ShortPtSpectrumBeforeCuts"), v0.pt());
              // Implementing best kzero cuts
              if (v0.v0cosPA() > kaonshSettingcosPA && v0.dcaV0daughters() < kaonshSettingdcav0dau && v0.v0radius() > kaonshSettingradius && std::abs(v0.dcapostopv()) > kaonshSettingdcapostopv && std::abs(v0.dcanegtopv()) > kaonshSettingdcanegtopv) {
                rPtAnalysis.fill(HIST("hMassK0ShortAllAfterCuts"), v0.mK0Short());
                rPtAnalysis.fill(HIST("hK0ShortReconstructedPtSpectrum"), v0.pt());
                for (int i = 0; i < 20; i++) {
                  // getting the pt value in #_# for and converting it to a number #.# for use, we get two values which correspond to the range of each bin
                  std::string pt1 = pthistos::kaonptbins[i];                // getting the lower string-value of the bin
                  std::string pt2 = pthistos::kaonptbins[i + 1];            // getting the higher string-value of the bin
                  size_t pos1 = pt1.find("_");                              // finding the "_" character of the lower string-value
                  size_t pos2 = pt2.find("_");                              // finding the "_" character of the higher string-value
                  pt1[pos1] = '.';                                          // changing the "_" character of the lower string-value to a "."
                  pt2[pos2] = '.';                                          // changing the "_" character of the higher string-value to a "."
                  const float ptlowervalue = std::stod(pt1);                // converting the lower string value to a double
                  const float pthighervalue = std::stod(pt2);               // converting the higher string value to a double
                  if (ptlowervalue <= v0.pt() && v0.pt() < pthighervalue) { // finding v0s with pt withing the range of our lower and higher value
                    pthistos::kaonPt[i]->Fill(v0.mK0Short());               // filling the 20 kaon namespace histograms
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
              // Implementing best lambda cuts
              if (v0.v0cosPA() > lambdaSettingcosPA && v0.dcaV0daughters() < lambdaSettingdcav0dau && v0.v0radius() > lambdaSettingradius && std::abs(v0.dcapostopv()) > lambdaSettingdcapostopv && std::abs(v0.dcanegtopv()) > lambdaSettingdcanegtopv) {
                rPtAnalysis.fill(HIST("hMassLambdaAllAfterCuts"), v0.mLambda());
                rPtAnalysis.fill(HIST("hLambdaReconstructedPtSpectrum"), v0.pt());
                for (int i = 0; i < 20; i++) {
                  // same as above with kzerosh we fill the 20 lambda namespace histograms within their Pt range
                  std::string pt1 = pthistos::lambdaPtBins[i];
                  std::string pt2 = pthistos::lambdaPtBins[i + 1];
                  size_t pos1 = pt1.find("_");
                  size_t pos2 = pt2.find("_");
                  pt1[pos1] = '.';
                  pt2[pos2] = '.';
                  const float ptlowervalue = std::stod(pt1);
                  const float pthighervalue = std::stod(pt2);
                  if (ptlowervalue <= v0.pt() && v0.pt() < pthighervalue) {
                    pthistos::lambdaPt[i]->Fill(v0.mLambda());
                  }
                }
              }
            }
          }
          // antilambda analysis
          if (antiLambdaAnalysis == true) {
            if (v0mcParticle.pdgCode() == -3122) { // antilambda matched
              rPtAnalysis.fill(HIST("hMassAntilambdaAll"), v0.mAntiLambda());
              rPtAnalysis.fill(HIST("hAntilambdaPtSpectrumBeforeCuts"), v0.pt());
              // Implementing best antilambda cuts
              if (v0.v0cosPA() > antilambdaSettingcosPA && v0.dcaV0daughters() < antilambdaSettingdcav0dau && v0.v0radius() > antilambdaSettingradius && std::abs(v0.dcapostopv()) > antilambdaSettingdcapostopv && std::abs(v0.dcanegtopv()) > antilambdaSettingdcanegtopv) {
                rPtAnalysis.fill(HIST("hMassAntilambdaAllAfterCuts"), v0.mAntiLambda());
                rPtAnalysis.fill(HIST("hAntilambdaReconstructedPtSpectrum"), v0.pt());
                for (int i = 0; i < 20; i++) {
                  // same as above with kzerosh and lambda we fill the 20 anti-lambda namespace histograms within their Pt range
                  std::string pt1 = pthistos::antiLambdaPtBins[i];
                  std::string pt2 = pthistos::antiLambdaPtBins[i + 1];
                  size_t pos1 = pt1.find("_");
                  size_t pos2 = pt2.find("_");
                  pt1[pos1] = '.';
                  pt2[pos2] = '.';
                  const float ptlowervalue = std::stod(pt1);
                  const float pthighervalue = std::stod(pt2);
                  if (ptlowervalue <= v0.pt() && v0.pt() < pthighervalue) {
                    pthistos::antiLambdaPt[i]->Fill(v0.mAntiLambda());
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
  void dataProcess(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const&,
                   aod::V0Datas const& V0s)
  {
    for (const auto& v0 : V0s) {
      // kzero analysis
      if (kzeroAnalysis == true) {
        // Filling the five Kzero invariant mass plots for different cuts (which are taken from namespace), for full explanation see the first kzero cut filling in the MC process
        rPtAnalysis.fill(HIST("hMassK0ShortAll"), v0.mK0Short());
        // Implementing best kzero cuts
        if (v0.v0cosPA() > kaonshSettingcosPA && v0.dcaV0daughters() < kaonshSettingdcav0dau && v0.v0radius() > kaonshSettingradius && std::abs(v0.dcapostopv()) > kaonshSettingdcapostopv && std::abs(v0.dcanegtopv()) > kaonshSettingdcanegtopv) {
          rPtAnalysis.fill(HIST("hMassK0ShortAllAfterCuts"), v0.mK0Short());
          for (int i = 0; i < 20; i++) { // same as above MC-process we fill the namespace histos with the kaon invariant mass of the particle within the pt range of the histo
            std::string pt1 = pthistos::kaonptbins[i];
            std::string pt2 = pthistos::kaonptbins[i + 1];
            size_t pos1 = pt1.find("_");
            size_t pos2 = pt2.find("_");
            pt1[pos1] = '.';
            pt2[pos2] = '.';
            const float ptlowervalue = std::stod(pt1);
            const float pthighervalue = std::stod(pt2);
            if (ptlowervalue < v0.pt() && v0.pt() < pthighervalue) {
              pthistos::kaonPt[i]->Fill(v0.mK0Short());
            }
          }
        }
      }
      // lambda analysis
      if (lambdaAnalysis == true) {
        // Filling the five lambda invariant mass plots for different cuts (which are taken from namespace), for full explanation see the first kzero cut filling in the MC process
        rPtAnalysis.fill(HIST("hMassLambdaAll"), v0.mLambda());
        // Implementing best lambda cuts
        if (v0.v0cosPA() > lambdaSettingcosPA && v0.dcaV0daughters() < lambdaSettingdcav0dau && v0.v0radius() > lambdaSettingradius && std::abs(v0.dcapostopv()) > lambdaSettingdcapostopv && std::abs(v0.dcanegtopv()) > lambdaSettingdcanegtopv) {
          rPtAnalysis.fill(HIST("hMassLambdaAllAfterCuts"), v0.mLambda());
          for (int i = 0; i < 20; i++) { // same as above MC-process we fill the namespace histos with the lambda invariant mass of the particle within the pt range of the histo
            std::string pt1 = pthistos::lambdaPtBins[i];
            std::string pt2 = pthistos::lambdaPtBins[i + 1];
            size_t pos1 = pt1.find("_");
            size_t pos2 = pt2.find("_");
            pt1[pos1] = '.';
            pt2[pos2] = '.';
            const float ptlowervalue = std::stod(pt1);
            const float pthighervalue = std::stod(pt2);
            if (ptlowervalue < v0.pt() && v0.pt() < pthighervalue) {
              pthistos::lambdaPt[i]->Fill(v0.mLambda());
            }
          }
        }
      }
      // anti-lambda analysis
      if (antiLambdaAnalysis == true) {
        // Filling the five Antilambda invariant mass plots for different cuts (which are taken from namespace), for full explanation see the first kzero cut filling in the MC process
        rPtAnalysis.fill(HIST("hMassAntilambdaAll"), v0.mAntiLambda());
        // implementing best antilambda cuts
        if (v0.v0cosPA() > antilambdaSettingcosPA && v0.dcaV0daughters() < antilambdaSettingdcav0dau && v0.v0radius() > antilambdaSettingradius && std::abs(v0.dcapostopv()) > antilambdaSettingdcapostopv && std::abs(v0.dcanegtopv()) > antilambdaSettingdcanegtopv) {
          rPtAnalysis.fill(HIST("hMassAntilambdaAllAfterCuts"), v0.mAntiLambda());
          for (int i = 0; i < 20; i++) { // same as above MC-process we fill the namespace histos with the antilambda invariant mass of the particle within the pt range of the histo
            std::string pt1 = pthistos::antiLambdaPtBins[i];
            std::string pt2 = pthistos::antiLambdaPtBins[i + 1];
            size_t pos1 = pt1.find("_");
            size_t pos2 = pt2.find("_");
            pt1[pos1] = '.';
            pt2[pos2] = '.';
            const float ptlowervalue = std::stod(pt1);
            const float pthighervalue = std::stod(pt2);
            if (ptlowervalue < v0.pt() && v0.pt() < pthighervalue) {
              pthistos::antiLambdaPt[i]->Fill(v0.mAntiLambda());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(V0PtInvMassPlots, genMCProcess, "Process Run 3 MC Generated", false);
  PROCESS_SWITCH(V0PtInvMassPlots, recMCProcess, "Process Run 3 MC", false);
  PROCESS_SWITCH(V0PtInvMassPlots, dataProcess, "Process Run 3 Data,", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<V0PtInvMassPlots>(cfgc)};
}
