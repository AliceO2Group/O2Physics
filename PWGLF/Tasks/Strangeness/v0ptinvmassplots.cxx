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
/// \brief V0 task for production of invariant mass plots for Optimised Topological Cuts Analysis
/// \author Nikolaos Karatzenis (nikolaos.karatzenis@cern.ch)
/// \author Roman Lietava (roman.lietava@cern.ch)

/*Description
This task creates 20 histograms that are filled with the V0 invariant mass under the K0, Lambda and Antilambda mass assumption
for different pt ranges (constituting bins), so 3x20=60 plots.The values are inserted as configurable strings for convinience.
Plots of the invariant masses at different stages of the analysis (ex. before and after the V0 cuts are enforced) and some pt distributions.
This analysis includes two processes, one for Real Data and one for MC Data switchable at the end of the code, only run one at a time*/

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"

// namespace to be used for pt plots and bins
namespace pthistos
{
std::shared_ptr<TH1> KaonPt[20];
static std::vector<std::string> kaonptbins;
std::shared_ptr<TH1> LambdaPt[20];
static std::vector<std::string> lambdaptbins;
std::shared_ptr<TH1> AntilambdaPt[20];
static std::vector<std::string> antilambdaptbins;
} // namespace pthistos
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct v0ptinvmassplots {
  // Histogram Registries
  HistogramRegistry rPtAnalysis{"PtAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKaonshMassPlots_per_PtBin{"KaonshMassPlots_per_PtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambdaMassPlots_per_PtBin{"LambdaMassPlots_per_PtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntilambdaMassPlots_per_PtBin{"AntilambdaMassPlots_per_PtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  Configurable<int> xaxisgenbins{"xaxisgenbins", 20, "Number of bins for Generated Pt Spectrum"};
  Configurable<float> xaxismingenbin{"xaxismingenbin", 0.0, "Minimum bin value of the Generated Pt Spectrum Plot"};
  Configurable<float> xaxismaxgenbin{"xaxismaxgenbin", 3.0, "Maximum bin value of the Generated Pt Spectrum Plot"};

  // Configurable Kaonsh Topological Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> kaonshsetting_dcav0dau{"kaonshsetting_dcav0dau", 100.0, "DCA V0 Daughters"};
  Configurable<float> kaonshsetting_dcapostopv{"kaonshsetting_dcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<float> kaonshsetting_dcanegtopv{"kaonshsetting_dcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<double> kaonshsetting_cospa{"kaonshsetting_cospa", 0.50, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> kaonshsetting_radius{"kaonshsetting_radius", 0.50, "v0radius"};

  // Configurable Lambda Topological Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> lambdasetting_dcav0dau{"lambdasetting_dcav0dau", 100.0, "DCA V0 Daughters"};
  Configurable<float> lambdasetting_dcapostopv{"lambdasetting_dcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<float> lambdasetting_dcanegtopv{"lambdasetting_dcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<double> lambdasetting_cospa{"lambdasetting_cospa", 0.50, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> lambdasetting_radius{"lambdasetting_radius", 0.50, "v0radius"};

  // Configurable Antilambda Topological Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> antilambdasetting_dcav0dau{"antilambdasetting_dcav0dau", 100.0, "DCA V0 Daughters"};
  Configurable<float> antilambdasetting_dcapostopv{"antilambdasetting_dcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<float> antilambdasetting_dcanegtopv{"antilambdasetting_dcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<double> antilambdasetting_cospa{"antilambdasetting_cospa", 0.50, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> antilambdasetting_radius{"antilambdasetting_radius", 0.50, "v0radius"};

  // Configurables for Specific V0s analysis
  Configurable<bool> kzerosh_analysis{"kzerosh_analysis", true, "Enable Kzerosh Pt Analysis"};
  Configurable<bool> lambda_analysis{"lambda_analysis", true, "Enable Lambda Pt Analysis"};
  Configurable<bool> antilambda_analysis{"antilambda_analysis", true, "Enable Antilambda Pt Analysis"};

  // Configurable string for Different Pt Bins
  Configurable<std::string> kzeroshsetting_pt_string{"kzerosetting_ptbins", {"0_0,0_15,0_3,0_45,0_6,0_75,0_9,1_05,1_2,1_35,1_5,1_65,1_8,1_95,2_1,2_25,2_4,2_55,2_7,2_85,3_0"}, "Kzero Pt Bin Values"};
  Configurable<std::string> lambdasetting_pt_string{"lambdasetting_ptbins", {"0_0,0_15,0_3,0_45,0_6,0_75,0_9,1_05,1_2,1_35,1_5,1_65,1_8,1_95,2_1,2_25,2_4,2_55,2_7,2_85,3_0"}, "Lambda Pt Bin Values"};
  Configurable<std::string> antilambdasetting_pt_string{"antilambdasetting_ptbins", {"0_0,0_15,0_3,0_45,0_6,0_75,0_9,1_05,1_2,1_35,1_5,1_65,1_8,1_95,2_1,2_25,2_4,2_55,2_7,2_85,3_0"}, "Antilambda Pt Bin Values"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec K0ShortMassAxis = {nBins, 0.45f, 0.55f, "#it{M} #pi^{+}#pi^{-} [GeV/#it{c}^{2}]"};
    AxisSpec LambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M} p^{+}#pi^{-} [GeV/#it{c}^{2}]"};
    AxisSpec AntiLambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M} p^{-}#pi^{+} [GeV/#it{c}^{2}]"};
    AxisSpec ptAxis = {nBins, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec GenptAxis = {xaxisgenbins, xaxismingenbin, xaxismaxgenbin, "#it{p}_{T} (GeV/#it{c})"};

    rPtAnalysis.add("hV0PtAll", "hV0PtAll", {HistType::kTH1F, {{nBins, 0.0f, 10.0f}}});

    // setting strings from configurable strings in order to manipulate them
    size_t commapos = 0;
    std::string token1;
    // Adding Kzerosh Histograms to registry
    if (kzerosh_analysis == true) {
      // getting the  bin values for the names of the plots for the five topological cuts
      std::string kzeroshsetting_ptbins = kzeroshsetting_pt_string;
      for (int i = 0; i < 21; i++) {                        // we have 21 pt values (for the 20 histos) as they are the ranges
        commapos = kzeroshsetting_ptbins.find(",");         // find comma that separates the values in the string
        token1 = kzeroshsetting_ptbins.substr(0, commapos); // store the substring (first individual value)
        pthistos::kaonptbins.push_back(token1);             //  fill the namespace with the value
        kzeroshsetting_ptbins.erase(0, commapos + 1);       // erase the value from the set string so it moves to the next
      }
      rPtAnalysis.add("hK0ShortReconstructedPtSpectrum", "hK0ShortReconstructedPtSpectrum", {HistType::kTH1F, {ptAxis}});
      rPtAnalysis.add("hMassK0ShortAll", "hMassK0ShortAll", {HistType::kTH1F, {K0ShortMassAxis}});
      rPtAnalysis.add("hK0ShortPtSpectrumBeforeCuts", "hK0ShortPtSpectrumBeforeCuts", {HistType::kTH1F, {ptAxis}});
      rPtAnalysis.add("hMassK0ShortAllAfterCuts", "hMassK0ShortAllAfterCuts", {HistType::kTH1F, {K0ShortMassAxis}});
      rPtAnalysis.add("hK0ShGeneratedPtSpectrum", "hK0ShGeneratedPtSpectrum", {HistType::kTH1F, {GenptAxis}});
      rPtAnalysis.add("hLambdaGeneratedPtSpectrum", "hLambdaGeneratedPtSpectrum", {HistType::kTH1F, {GenptAxis}});
      rPtAnalysis.add("hAntilambdaGeneratedPtSpectrum", "hAntilambdaGeneratedPtSpectrum", {HistType::kTH1F, {GenptAxis}});
      for (int i = 0; i < 20; i++) {
        pthistos::KaonPt[i] = rKaonshMassPlots_per_PtBin.add<TH1>(fmt::format("hPt_from_{0}_to_{1}", pthistos::kaonptbins[i], pthistos::kaonptbins[i + 1]).data(), fmt::format("hPt from {0} to {1}", pthistos::kaonptbins[i], pthistos::kaonptbins[i + 1]).data(), {HistType::kTH1D, {{K0ShortMassAxis}}});
      }
    }
    // Adding Lambda Histograms
    if (lambda_analysis == true) {
      // same method as in Kzerosh above
      std::string lambdasetting_ptbins = lambdasetting_pt_string;
      for (int i = 0; i < 21; i++) {
        commapos = lambdasetting_ptbins.find(",");
        token1 = lambdasetting_ptbins.substr(0, commapos);
        pthistos::lambdaptbins.push_back(token1);
        lambdasetting_ptbins.erase(0, commapos + 1);
      }
      rPtAnalysis.add("hLambdaReconstructedPtSpectrum", "hLambdaReconstructedPtSpectrum", {HistType::kTH1F, {ptAxis}});
      rPtAnalysis.add("hMassLambdaAll", "hMassLambdaAll", {HistType::kTH1F, {LambdaMassAxis}});
      rPtAnalysis.add("hLambdaPtSpectrumBeforeCuts", "hLambdaPtSpectrumBeforeCuts", {HistType::kTH1F, {ptAxis}});
      rPtAnalysis.add("hMassLambdaAllAfterCuts", "hMassLambdaAllAfterCuts", {HistType::kTH1F, {LambdaMassAxis}});
      for (int i = 0; i < 20; i++) {
        pthistos::LambdaPt[i] = rLambdaMassPlots_per_PtBin.add<TH1>(fmt::format("hPt_from_{0}_to_{1}", pthistos::lambdaptbins[i], pthistos::lambdaptbins[i + 1]).data(), fmt::format("hPt from {0} to {1}", pthistos::lambdaptbins[i], pthistos::lambdaptbins[i + 1]).data(), {HistType::kTH1D, {{LambdaMassAxis}}});
      }
    }
    // Adding Antilambda Histograms
    if (antilambda_analysis == true) {
      // same method as in Lambda and Kzerosh above
      std::string antilambdasetting_ptbins = antilambdasetting_pt_string;
      for (int i = 0; i < 21; i++) {
        commapos = antilambdasetting_ptbins.find(",");
        token1 = antilambdasetting_ptbins.substr(0, commapos);
        pthistos::antilambdaptbins.push_back(token1);
        antilambdasetting_ptbins.erase(0, commapos + 1);
      }
      rPtAnalysis.add("hAntilambdaReconstructedPtSpectrum", "hAntilambdaReconstructedPtSpectrum", {HistType::kTH1F, {ptAxis}});
      rPtAnalysis.add("hMassAntilambdaAll", "hMassAntilambdaAll", {HistType::kTH1F, {AntiLambdaMassAxis}});
      rPtAnalysis.add("hAntilambdaPtSpectrumBeforeCuts", "hAntilambdaPtSpectrumBeforeCuts", {HistType::kTH1F, {ptAxis}});
      rPtAnalysis.add("hMassAntilambdaAllAfterCuts", "hMassAntilambdaAllAfterCuts", {HistType::kTH1F, {AntiLambdaMassAxis}});
      for (int i = 0; i < 20; i++) {
        pthistos::AntilambdaPt[i] = rAntilambdaMassPlots_per_PtBin.add<TH1>(fmt::format("hPt_from_{0}_to_{1}", pthistos::antilambdaptbins[i], pthistos::antilambdaptbins[i + 1]).data(), fmt::format("hPt from {0} to {1}", pthistos::antilambdaptbins[i], pthistos::antilambdaptbins[i + 1]).data(), {HistType::kTH1D, {{AntiLambdaMassAxis}}});
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
  void GenMCprocess(soa::Filtered<aod::McCollisions>::iterator const&,
                    const soa::SmallGroups<soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>>&,
                    aod::McParticles const& mcParticles)
  {
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary()) {
        if (TMath::Abs(mcParticle.y()) < 0.5f) {
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
  void RecMCprocess(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const&,
                    soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s,
                    DaughterTracks const&, // no need to define a variable for tracks, if we don't access them directly
                    aod::McParticles const& mcParticles)
  {
    for (const auto& v0 : V0s) {
      rPtAnalysis.fill(HIST("hV0PtAll"), v0.pt());
      // Checking that the V0 is a true K0s/Lambdas/Antilambdas and then filling the parameter histograms and the invariant mass plots for different cuts (which are taken from namespace)
      if (v0.has_mcParticle()) {
        auto v0mcParticle = v0.mcParticle();
        if (v0mcParticle.isPhysicalPrimary()) {
          if (kzerosh_analysis == true) {
            if (v0mcParticle.pdgCode() == 310) { // kzero matched
              rPtAnalysis.fill(HIST("hMassK0ShortAll"), v0.mK0Short());
              rPtAnalysis.fill(HIST("hK0ShortPtSpectrumBeforeCuts"), v0.pt());
              // Implementing best kzero cuts
              if (v0.v0cosPA() > kaonshsetting_cospa && v0.dcaV0daughters() < kaonshsetting_dcav0dau && v0.v0radius() > kaonshsetting_radius && TMath::Abs(v0.dcapostopv()) > kaonshsetting_dcapostopv && TMath::Abs(v0.dcanegtopv()) > kaonshsetting_dcanegtopv) {
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
                    pthistos::KaonPt[i]->Fill(v0.mK0Short());               // filling the 20 kaon namespace histograms
                  }
                }
              }
            }
          }
          // lambda analysis
          if (lambda_analysis == true) {
            if (v0mcParticle.pdgCode() == 3122) { // lambda matched
              rPtAnalysis.fill(HIST("hMassLambdaAll"), v0.mLambda());
              rPtAnalysis.fill(HIST("hLambdaPtSpectrumBeforeCuts"), v0.pt());
              // Implementing best lambda cuts
              if (v0.v0cosPA() > lambdasetting_cospa && v0.dcaV0daughters() < lambdasetting_dcav0dau && v0.v0radius() > lambdasetting_radius && TMath::Abs(v0.dcapostopv()) > lambdasetting_dcapostopv && TMath::Abs(v0.dcanegtopv()) > lambdasetting_dcanegtopv) {
                rPtAnalysis.fill(HIST("hMassLambdaAllAfterCuts"), v0.mLambda());
                rPtAnalysis.fill(HIST("hLambdaReconstructedPtSpectrum"), v0.pt());
                for (int i = 0; i < 20; i++) {
                  // same as above with kzerosh we fill the 20 lambda namespace histograms within their Pt range
                  std::string pt1 = pthistos::lambdaptbins[i];
                  std::string pt2 = pthistos::lambdaptbins[i + 1];
                  size_t pos1 = pt1.find("_");
                  size_t pos2 = pt2.find("_");
                  pt1[pos1] = '.';
                  pt2[pos2] = '.';
                  const float ptlowervalue = std::stod(pt1);
                  const float pthighervalue = std::stod(pt2);
                  if (ptlowervalue <= v0.pt() && v0.pt() < pthighervalue) {
                    pthistos::LambdaPt[i]->Fill(v0.mLambda());
                  }
                }
              }
            }
          }
          // antilambda analysis
          if (antilambda_analysis == true) {
            if (v0mcParticle.pdgCode() == -3122) { // antilambda matched
              rPtAnalysis.fill(HIST("hMassAntilambdaAll"), v0.mAntiLambda());
              rPtAnalysis.fill(HIST("hAntilambdaPtSpectrumBeforeCuts"), v0.pt());
              // Implementing best antilambda cuts
              if (v0.v0cosPA() > antilambdasetting_cospa && v0.dcaV0daughters() < antilambdasetting_dcav0dau && v0.v0radius() > antilambdasetting_radius && TMath::Abs(v0.dcapostopv()) > antilambdasetting_dcapostopv && TMath::Abs(v0.dcanegtopv()) > antilambdasetting_dcanegtopv) {
                rPtAnalysis.fill(HIST("hMassAntilambdaAllAfterCuts"), v0.mAntiLambda());
                rPtAnalysis.fill(HIST("hAntilambdaReconstructedPtSpectrum"), v0.pt());
                for (int i = 0; i < 20; i++) {
                  // same as above with kzerosh and lambda we fill the 20 anti-lambda namespace histograms within their Pt range
                  std::string pt1 = pthistos::antilambdaptbins[i];
                  std::string pt2 = pthistos::antilambdaptbins[i + 1];
                  size_t pos1 = pt1.find("_");
                  size_t pos2 = pt2.find("_");
                  pt1[pos1] = '.';
                  pt2[pos2] = '.';
                  const float ptlowervalue = std::stod(pt1);
                  const float pthighervalue = std::stod(pt2);
                  if (ptlowervalue <= v0.pt() && v0.pt() < pthighervalue) {
                    pthistos::AntilambdaPt[i]->Fill(v0.mAntiLambda());
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
  void Dataprocess(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const&,
                   aod::V0Datas const& V0s)
  {
    for (const auto& v0 : V0s) {
      // kzero analysis
      if (kzerosh_analysis == true) {
        // Filling the five Kzero invariant mass plots for different cuts (which are taken from namespace), for full explanation see the first kzero cut filling in the MC process
        rPtAnalysis.fill(HIST("hMassK0ShortAll"), v0.mK0Short());
        // Implementing best kzero cuts
        if (v0.v0cosPA() > kaonshsetting_cospa && v0.dcaV0daughters() < kaonshsetting_dcav0dau && v0.v0radius() > kaonshsetting_radius && TMath::Abs(v0.dcapostopv()) > kaonshsetting_dcapostopv && TMath::Abs(v0.dcanegtopv()) > kaonshsetting_dcanegtopv) {
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
              pthistos::KaonPt[i]->Fill(v0.mK0Short());
            }
          }
        }
      }
      // lambda analysis
      if (lambda_analysis == true) {
        // Filling the five lambda invariant mass plots for different cuts (which are taken from namespace), for full explanation see the first kzero cut filling in the MC process
        rPtAnalysis.fill(HIST("hMassLambdaAll"), v0.mLambda());
        // Implementing best lambda cuts
        if (v0.v0cosPA() > lambdasetting_cospa && v0.dcaV0daughters() < lambdasetting_dcav0dau && v0.v0radius() > lambdasetting_radius && TMath::Abs(v0.dcapostopv()) > lambdasetting_dcapostopv && TMath::Abs(v0.dcanegtopv()) > lambdasetting_dcanegtopv) {
          rPtAnalysis.fill(HIST("hMassLambdaAllAfterCuts"), v0.mLambda());
          for (int i = 0; i < 20; i++) { // same as above MC-process we fill the namespace histos with the lambda invariant mass of the particle within the pt range of the histo
            std::string pt1 = pthistos::lambdaptbins[i];
            std::string pt2 = pthistos::lambdaptbins[i + 1];
            size_t pos1 = pt1.find("_");
            size_t pos2 = pt2.find("_");
            pt1[pos1] = '.';
            pt2[pos2] = '.';
            const float ptlowervalue = std::stod(pt1);
            const float pthighervalue = std::stod(pt2);
            if (ptlowervalue < v0.pt() && v0.pt() < pthighervalue) {
              pthistos::LambdaPt[i]->Fill(v0.mLambda());
            }
          }
        }
      }
      // anti-lambda analysis
      if (antilambda_analysis == true) {
        // Filling the five Antilambda invariant mass plots for different cuts (which are taken from namespace), for full explanation see the first kzero cut filling in the MC process
        rPtAnalysis.fill(HIST("hMassAntilambdaAll"), v0.mAntiLambda());
        // implementing best antilambda cuts
        if (v0.v0cosPA() > antilambdasetting_cospa && v0.dcaV0daughters() < antilambdasetting_dcav0dau && v0.v0radius() > antilambdasetting_radius && TMath::Abs(v0.dcapostopv()) > antilambdasetting_dcapostopv && TMath::Abs(v0.dcanegtopv()) > antilambdasetting_dcanegtopv) {
          rPtAnalysis.fill(HIST("hMassAntilambdaAllAfterCuts"), v0.mAntiLambda());
          for (int i = 0; i < 20; i++) { // same as above MC-process we fill the namespace histos with the antilambda invariant mass of the particle within the pt range of the histo
            std::string pt1 = pthistos::antilambdaptbins[i];
            std::string pt2 = pthistos::antilambdaptbins[i + 1];
            size_t pos1 = pt1.find("_");
            size_t pos2 = pt2.find("_");
            pt1[pos1] = '.';
            pt2[pos2] = '.';
            const float ptlowervalue = std::stod(pt1);
            const float pthighervalue = std::stod(pt2);
            if (ptlowervalue < v0.pt() && v0.pt() < pthighervalue) {
              pthistos::AntilambdaPt[i]->Fill(v0.mAntiLambda());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(v0ptinvmassplots, GenMCprocess, "Process Run 3 MC Generated", false);
  PROCESS_SWITCH(v0ptinvmassplots, RecMCprocess, "Process Run 3 MC", false);
  PROCESS_SWITCH(v0ptinvmassplots, Dataprocess, "Process Run 3 Data,", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0ptinvmassplots>(cfgc)};
}
