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

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include "CommonConstants/PhysicsConstants.h"
#include "CommonUtils/StringUtils.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "TPDGCode.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

// namespace to be used for pt plots and bins
namespace pthistos
{
std::vector<std::shared_ptr<TH2>> kaonPt;
static std::vector<std::string> kaonPtBins;
std::vector<std::shared_ptr<TH2>> lambdaPt;
static std::vector<std::string> lambdaPtBins;
std::vector<std::shared_ptr<TH2>> antilambdaPt;
static std::vector<std::string> antilambdaPtBins;
std::vector<std::shared_ptr<TH2>> kaonSplit;
std::vector<std::shared_ptr<TH2>> lambdaSplit;
std::vector<std::shared_ptr<TH2>> antilambdaSplit;
} // namespace pthistos
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct V0PtInvMassPlots {
  // Histogram Registries
  HistogramRegistry rPtAnalysis{"PtAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rK0shMassPlotsPerPtBin{"K0shMassPlotsPerPtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambdaMassPlotsPerPtBin{"LambdaMassPlotsPerPtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntilambdaMassPlotsPerPtBin{"AntiLambdaMassPlotsPerPtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rK0shSplitMassPlotsPerPtBin{"K0shSplitMassPlotsPerPtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambdaSplitMassPlotsPerPtBin{"LambdaSplitMassPlotsPerPtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntilambdaSplitMassPlotsPerPtBin{"AntiLambdaSplitMassPlotsPerPtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rFeeddownMatrices{"FeeddownMatrices", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rMCCorrections{"MCCorrections", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  Configurable<int> nBinsArmenteros{"nBinsArmenteros", 500, "N bins in Armenteros histos"};

  // Configurables for Cuts
  Configurable<float> cutZVertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<float> nSigmaTPCPion{"nSigmaTPCPion", 4, "nSigmaTPCPion"};
  Configurable<float> nSigmaTPCProton{"nSigmaTPCProton", 4, "nSigmaTPCProton"};
  Configurable<float> compv0masscut{"compv0masscut", 0.01, "CompetitiveV0masscut (GeV)"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};
  Configurable<float> rapidityCut{"rapidityCut", 0.5, "V0 Rapidity Window"};
  Configurable<float> itsMinHits{"itsMinHits", 1.0, "Minimum Hits of Daughter Tracks in the ITS"};

  // Configurables switches for event selection
  Configurable<bool> dosel8{"dosel8", true, "Enable sel8 event selection"};
  Configurable<bool> doNoTimeFrameBorder{"doNoTimeFrameBorder", true, "Enable NoTimeFrameBorder event selection"};
  Configurable<bool> doNoITSROFrameBorder{"doNoITSROFrameBorder", true, "Enable NoITSROFrameBorder event selection"};
  Configurable<bool> doIsTriggerTVX{"doIsTriggerTVX", true, "Enable IsTriggerTVX event selection"};
  Configurable<bool> docutZVertex{"docutZVertex", true, "Enable cutZVertex event selection"};
  Configurable<bool> doIsVertexTOFmatched{"doIsVertexTOFmatched", true, "Enable IsVertexTOFmatched event selection"};
  Configurable<bool> doNoSameBunchPileup{"doNoSameBunchPileup", true, "Enable NoSameBunchPileup event selection"};
  Configurable<bool> doIsVertexITSTPC{"doIsVertexITSTPC", true, "Enable IsVertexITSTPC event selection"};
  Configurable<bool> doisInelGt0{"doisInelGt0", true, "Enable isInelGt0 event selection"};

  // Configurables switches for v0 selection
  Configurable<bool> doRapidityCut{"doRapidityCut", true, "Enable rapidity v0 selection"};
  Configurable<bool> doDaughterPseudorapidityCut{"doDaughterPseudorapidityCut", true, "Enable Daughter pseudorapidity v0 selection"};
  Configurable<bool> doisNotITSAfterburner{"doisNotITSAfterburner", true, "Enable Tracks do not come from Afterburner"};
  Configurable<bool> doitsMinHits{"doitsMinHits", true, "Enable ITS Minimum hits"};

  // Configurables switches for K0sh selection
  Configurable<bool> dotruthk0sh{"dotruthk0sh", true, "Enable K0sh MC Matching"};
  Configurable<bool> doK0shTPCPID{"doK0shTPCPID", true, "Enable K0sh TPC PID"};
  Configurable<bool> doK0shcomptmasscut{"doK0shcomptmasscut", true, "Enable K0sh Competitive V0 Mass Cut"};
  Configurable<bool> doK0shMaxct{"doK0shMaxct", true, "Enable K0sh Max ct Cut"};
  Configurable<bool> doK0shArmenterosCut{"doK0shArmenterosCut", true, "Enable K0sh Armenteros Cut"};
  Configurable<bool> doK0shcosPACut{"doK0shcosPACut", true, "Enable K0sh cosPA Topological Cut"};
  Configurable<bool> doK0shDCAdauCut{"doK0shDCAdauCut", true, "Enable K0sh DCA daughters Topological Cut"};
  Configurable<bool> doK0shv0radiusCut{"doK0shv0radiusCut", true, "Enable K0sh v0radius Topological Cut"};
  Configurable<bool> doK0shdcaposdautopv{"doK0shdcaposdautopv", true, "Enable K0sh DCA pos daughter to PV Topological Cut"};
  Configurable<bool> doK0shdcanegdautopv{"doK0shdcanegdautopv", true, "Enable K0sh DCA neg daughter to PV Topological Cut"};

  // Configurables switches for Lambda selection
  Configurable<bool> dotruthLambda{"dotruthLambda", true, "Enable Lambda MC Matching"};
  Configurable<bool> doLambdaTPCPID{"doLambdaTPCPID", true, "Enable Lambda TPC PID"};
  Configurable<bool> doLambdacomptmasscut{"doLambdacomptmasscut", true, "Enable Lambda Competitive V0 Mass Cut"};
  Configurable<bool> doLambdaMaxct{"doLambdaMaxct", true, "Enable Lambda Max ct Cut"};
  Configurable<bool> doLambdaArmenterosCut{"doLambdaArmenterosCut", true, "Enable Lambda Armenteros Cut"};
  Configurable<bool> doLambdacosPACut{"doLambdacosPACut", true, "Enable Lambda cosPA Topological Cut"};
  Configurable<bool> doLambdaDCAdauCut{"doLambdaDCAdauCut", true, "Enable Lambda DCA daughters Topological Cut"};
  Configurable<bool> doLambdav0radiusCut{"doLambdav0radiusCut", true, "Enable Lambda v0radius Topological Cut"};
  Configurable<bool> doLambdadcaposdautopv{"doLambdadcaposdautopv", true, "Enable Lambda DCA pos daughter to PV Topological Cut"};
  Configurable<bool> doLambdadcanegdautopv{"doLambdadcanegdautopv", true, "Enable Lambda DCA neg daughter to PV Topological Cut"};

  // Configurables switches for Lambda selection
  Configurable<bool> dotruthAntilambda{"dotruthAntilambda", true, "Enable Antilambda MC Matching"};
  Configurable<bool> doAntilambdaTPCPID{"doAntilambdaTPCPID", true, "Enable Antilambda TPC PID"};
  Configurable<bool> doAntilambdacomptmasscut{"doAntilambdacomptmasscut", true, "Enable Antilambda Competitive V0 Mass Cut"};
  Configurable<bool> doAntilambdaMaxct{"doAntilambdaMaxct", true, "Enable Antilambda Max ct Cut"};
  Configurable<bool> doAntilambdaArmenterosCut{"doAntilambdaArmenterosCut", true, "Enable Antilambda Armenteros Cut"};
  Configurable<bool> doAntilambdacosPACut{"doAntilambdacosPACut", true, "Enable Antilambda cosPA Topological Cut"};
  Configurable<bool> doAntilambdaDCAdauCut{"doAntilambdaDCAdauCut", true, "Enable Antilambda DCA daughters Topological Cut"};
  Configurable<bool> doAntilambdav0radiusCut{"doAntilambdav0radiusCut", true, "Enable Antilambda v0radius Topological Cut"};
  Configurable<bool> doAntilambdadcaposdautopv{"doAntilambdadcaposdautopv", true, "Enable Antilambda DCA pos daughter to PV Topological Cut"};
  Configurable<bool> doAntilambdadcanegdautopv{"doAntilambdadcanegdautopv", true, "Enable Antilambda DCA neg daughter to PV Topological Cut"};

  // Configurable K0sh Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> k0shSettingdcav0dau{"k0shSettingdcav0dau", 0.3, "DCA V0 Daughters"};
  Configurable<float> k0shSettingdcapostopv{"k0shSettingdcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<float> k0shSettingdcanegtopv{"k0shSettingdcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<double> k0shSettingcosPA{"k0shSettingcosPA", 0.98, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> k0shSettingradius{"k0shSettingradius", 0.50, "v0radius"};
  Configurable<float> k0shmaxct{"k0shmaxct", 20.00, "K0sh maximum ct value"};
  Configurable<float> k0shparamArmenterosCut{"k0shparamArmenterosCut", 0.2, "K0sh Armenteros Cut on parameter"};

  // Configurable Lambda Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> lambdaSettingdcav0dau{"lambdaSettingdcav0dau", 0.3, "DCA V0 Daughters"};
  Configurable<float> lambdaSettingdcapostopv{"lambdaSettingdcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<float> lambdaSettingdcanegtopv{"lambdaSettingdcanegtopv", 0.09, "DCA Neg To PV"};
  Configurable<double> lambdaSettingcosPA{"lambdaSettingcosPA", 0.98, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> lambdaSettingradius{"lambdaSettingradius", 0.50, "v0radius"};
  Configurable<float> lambdamaxct{"lambdamaxct", 30.00, "Lambda maximum ct value"};
  Configurable<float> lambdaparamArmenterosCut{"lambdaparamArmenterosCut", 0.2, "Lambda Armenteros Cut on parameter"};

  // Configurable Antilambda Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> antilambdaSettingdcav0dau{"antilambdaSettingdcav0dau", 0.3, "DCA V0 Daughters"};
  Configurable<float> antilambdaSettingdcapostopv{"antilambdaSettingdcapostopv", 0.09, "DCA Pos To PV"};
  Configurable<float> antilambdaSettingdcanegtopv{"antilambdaSettingdcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<double> antilambdaSettingcosPA{"antilambdaSettingcosPA", 0.98, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> antilambdaSettingradius{"antilambdaSettingradius", 0.50, "v0radius"};
  Configurable<float> antilambdamaxct{"antilambdamaxct", 30.00, "Antilambda maximum ct value"};
  Configurable<float> antilambdaparamArmenterosCut{"antilambdaparamArmenterosCut", 0.2, "Antilambda Armenteros Cut on parameter"};

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
    // tokenise strings into individual values
    pthistos::kaonPtBins = o2::utils::Str::tokenize(kzeroSettingPtBinsString, ',');
    pthistos::lambdaPtBins = o2::utils::Str::tokenize(lambdaSettingPtBinsString, ',');
    pthistos::antilambdaPtBins = o2::utils::Str::tokenize(antilambdaSettingPtBinsString, ',');

    // Calculate number of histograms for each particle type
    int nKaonHistograms = pthistos::kaonPtBins.size() - 1;
    int nLambdaHistograms = pthistos::lambdaPtBins.size() - 1;
    int nAntilambdaHistograms = pthistos::antilambdaPtBins.size() - 1;

    pthistos::kaonPt.resize(nKaonHistograms);                // number of Kaon Pt histograms to expect
    pthistos::lambdaPt.resize(nLambdaHistograms);            // number of Lambda histograms to expect
    pthistos::antilambdaPt.resize(nAntilambdaHistograms);    // number of Antilambda histograms to expect
    pthistos::kaonSplit.resize(nKaonHistograms);             // number of Kaon Split Pt histograms to expect
    pthistos::lambdaSplit.resize(nLambdaHistograms);         // number of Lambda Split Pt histograms to expect
    pthistos::antilambdaSplit.resize(nAntilambdaHistograms); // number of Antilambda Split Pt histograms to expect

    // initialize and convert tokenized strings into vector of doubles for AxisSpec
    std::vector<double> kaonptedgevalues(pthistos::kaonPtBins.size());
    std::vector<double> lambdaptedgevalues(pthistos::lambdaPtBins.size());
    std::vector<double> antilambdaptedgevalues(pthistos::antilambdaPtBins.size());
    for (size_t i = 0; i < pthistos::kaonPtBins.size(); i++) {
      kaonptedgevalues[i] = std::stod(pthistos::kaonPtBins[i]);
    }
    for (size_t i = 0; i < pthistos::lambdaPtBins.size(); i++) {
      lambdaptedgevalues[i] = std::stod(pthistos::lambdaPtBins[i]);
    }
    for (size_t i = 0; i < pthistos::antilambdaPtBins.size(); i++) {
      antilambdaptedgevalues[i] = std::stod(pthistos::antilambdaPtBins[i]);
    }

    // Axes
    AxisSpec k0ShortMassAxis = {nBins, 0.45f, 0.55f, "#it{M} #pi^{+}#pi^{-} [GeV/#it{c}^{2}]"};
    AxisSpec lambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M} p^{+}#pi^{-} [GeV/#it{c}^{2}]"};
    AxisSpec antiLambdaMassAxis = {nBins, 1.085f, 1.145f, "#it{M} p^{-}#pi^{+} [GeV/#it{c}^{2}]"};
    AxisSpec k0ShortPtAxis = {kaonptedgevalues, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec lambdaPtAxis = {lambdaptedgevalues, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec antilambdaPtAxis = {antilambdaptedgevalues, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec centAxis = {100, 0.0f, 100.0f, "#it{Centrality} (%)"};
    AxisSpec armenterosQtAxis = {nBinsArmenteros, 0.0f, 0.3f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec armenterosasymAxis = {nBinsArmenteros, -1.f, 1.f, "#it{p}^{+}_{||}-#it{p}^{-}_{||}/#it{p}^{+}_{||}+#it{p}^{-}_{||}"};
    AxisSpec vertexZAxis = {nBins, -11.0f, 11.0f, "vrtx_{Z} [cm]"};
    AxisSpec partCutsAxis{10, 0.0f, 10.0f, "Cut index"};

    std::vector<std::string> kaonhistvalue(nKaonHistograms + 1);
    std::vector<std::string> lambdahistvalue(nLambdaHistograms + 1);
    std::vector<std::string> antilambdahistvalue(nAntilambdaHistograms + 1);
    // K0short Histogram Pt Bin Edges (and Split)
    for (int i = 0; i < nKaonHistograms + 1; i++) {    // Histos won't accept "." character so converting it to "_"
      std::string kaonptbin = pthistos::kaonPtBins[i]; // getting the value of the bin edge
      size_t pos = kaonptbin.find(".");                // finding the "." character
      kaonptbin[pos] = '_';                            // changing the "." character of the string-value to a "_"
      kaonhistvalue[i] = kaonptbin;                    // filling bin edges list
    }
    // Lambda Histograms Pt Bin Edges (same as K0s above)
    for (int i = 0; i < nLambdaHistograms + 1; i++) {
      std::string lambdaptbin = pthistos::lambdaPtBins[i];
      size_t pos = lambdaptbin.find(".");
      lambdaptbin[pos] = '_';
      lambdahistvalue[i] = lambdaptbin;
    }
    // AntiLambda Histograms Pt Bin Edges (same as K0s above)
    for (int i = 0; i < nAntilambdaHistograms + 1; i++) {
      std::string antilambdaPtbin = pthistos::antilambdaPtBins[i];
      size_t pos = antilambdaPtbin.find(".");
      antilambdaPtbin[pos] = '_';
      antilambdahistvalue[i] = antilambdaPtbin;
    }

    // General Plots
    rPtAnalysis.add("hNEvents", "hNEvents", {HistType::kTH2D, {{7, 0.f, 7.f}, centAxis}});
    rPtAnalysis.add("hNRecEvents", "hNRecEvents", {HistType::kTH2D, {{1, 0.f, 1.f}, centAxis}});
    rPtAnalysis.add("hNV0s", "hNV0s", {HistType::kTH2D, {{10, 0.f, 10.f}, centAxis}});
    rPtAnalysis.add("hNK0sh", "hNK0sh", {HistType::kTH2D, {{11, 0.f, 11.f}, centAxis}});
    rPtAnalysis.add("hNLambda", "hNLambda", {HistType::kTH2D, {{11, 0.f, 11.f}, centAxis}});
    rPtAnalysis.add("hNAntiLambda", "hNAntiLambda", {HistType::kTH2D, {{11, 0.f, 11.f}, centAxis}});
    rPtAnalysis.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rPtAnalysis.add("hArmenterosPodolanskiPlot", "hArmenterosPodolanskiPlot", {HistType::kTH2F, {{armenterosasymAxis}, {armenterosQtAxis}}});
    rPtAnalysis.add("hV0EtaDaughters", "hV0EtaDaughters", {HistType::kTH1F, {{nBins, -1.2f, 1.2f}}});
    rPtAnalysis.add("V0Rapidity", "V0Rapidity", {HistType::kTH1F, {{nBins, -1.0f, 1.0f}}});

    // Adding Kzerosh Histograms to registry
    if (kzeroAnalysis == true) {
      rPtAnalysis.add("hMassK0ShortvsCuts", "hMassK0ShortvsCuts", {HistType::kTH2F, {{partCutsAxis}, {k0ShortMassAxis}}});
      rPtAnalysis.add("hArmenterosPodolanskiPlotK0sh", "hArmenterosPodolanskiPlotK0sh", {HistType::kTH2F, {{armenterosasymAxis}, {armenterosQtAxis}}});
      rPtAnalysis.add("hNSigmaPosPionFromK0s", "hNSigmaPosPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {k0ShortPtAxis}}});
      rPtAnalysis.add("hNSigmaNegPionFromK0s", "hNSigmaNegPionFromK0s", {HistType::kTH2F, {{100, -5.f, 5.f}, {k0ShortPtAxis}}});
      rPtAnalysis.add("hK0shV0radius", "hK0shV0radius", {HistType::kTH1F, {{nBins, 0.0f, 50.0f}}});
      rPtAnalysis.add("hK0shcosPA", "hK0shcosPA", {HistType::kTH1F, {{nBins, 0.95f, 1.0f}}});
      rPtAnalysis.add("hK0shDCAV0Daughters", "hK0shDCAV0Daughters", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      rPtAnalysis.add("hK0shDCAPosDaughter", "hK0shDCAPosDaughter", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      rPtAnalysis.add("hK0shDCANegDaughter", "hK0shDCANegDaughter", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      for (int i = 0; i < nKaonHistograms; i++) {
        pthistos::kaonPt[i] = rK0shMassPlotsPerPtBin.add<TH2>(fmt::format("hK0shPt_vs_Cent_from_{0}_to_{1}", kaonhistvalue[i], kaonhistvalue[i + 1]).c_str(), fmt::format("K0s mass vs centrality, pT from {0} to {1}", kaonhistvalue[i], kaonhistvalue[i + 1]).c_str(), {HistType::kTH2D, {{k0ShortMassAxis}, {centAxis}}});
        pthistos::kaonSplit[i] = rK0shSplitMassPlotsPerPtBin.add<TH2>(fmt::format("hK0shSplitPt_vs_Cent_from_{0}_to_{1}", kaonhistvalue[i], kaonhistvalue[i + 1]).c_str(), fmt::format("Split K0s mass vs centrality, pT from {0} to {1}", kaonhistvalue[i], kaonhistvalue[i + 1]).c_str(), {HistType::kTH2D, {{k0ShortMassAxis}, {centAxis}}});
      }
      rFeeddownMatrices.add("hK0shFeeddownMatrix", "hK0shFeeddownMatrix", {HistType::kTH3D, {k0ShortPtAxis, k0ShortPtAxis, centAxis}});
      rFeeddownMatrices.add("hK0shPhiFeeddownMatrix", "hK0shPhiFeeddownMatrix", {HistType::kTH3D, {k0ShortPtAxis, k0ShortPtAxis, centAxis}});
    }
    // Adding Lambda Histograms
    if (lambdaAnalysis == true) {
      // same method as in Kzerosh above
      rPtAnalysis.add("hMassLambdavsCuts", "hMassLambdavsCuts", {HistType::kTH2F, {{partCutsAxis}, {k0ShortMassAxis}}});
      rPtAnalysis.add("hArmenterosPodolanskiPlotLambda", "hArmenterosPodolanskiPlotLambda", {HistType::kTH2F, {{armenterosasymAxis}, {armenterosQtAxis}}});
      rPtAnalysis.add("hNSigmaPosProtonFromLambdas", "hNSigmaPosProtonFromLambdas", {HistType::kTH2F, {{100, -5.f, 5.f}, {lambdaPtAxis}}});
      rPtAnalysis.add("hNSigmaNegPionFromLambdas", "hNSigmaNegPionFromLambdas", {HistType::kTH2F, {{100, -5.f, 5.f}, {lambdaPtAxis}}});
      rPtAnalysis.add("hLambdaV0radius", "hLambdaV0radius", {HistType::kTH1F, {{nBins, 0.0f, 50.0f}}});
      rPtAnalysis.add("hLambdacosPA", "hLambdacosPA", {HistType::kTH1F, {{nBins, 0.95f, 1.0f}}});
      rPtAnalysis.add("hLambdaDCAV0Daughters", "hLambdaDCAV0Daughters", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      rPtAnalysis.add("hLambdaDCAPosDaughter", "hLambdaDCAPosDaughter", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      rPtAnalysis.add("hLambdaDCANegDaughter", "hLambdaDCANegDaughter", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});

      for (int i = 0; i < nLambdaHistograms; i++) {
        pthistos::lambdaPt[i] = rLambdaMassPlotsPerPtBin.add<TH2>(fmt::format("hLambdaPt_vs_Cent_from_{0}_to_{1}", lambdahistvalue[i], lambdahistvalue[i + 1]).c_str(), fmt::format("Lambda mass vs centrality, pT from {0} to {1}", lambdahistvalue[i], lambdahistvalue[i + 1]).c_str(), {HistType::kTH2D, {{lambdaMassAxis}, {centAxis}}});
        pthistos::lambdaSplit[i] = rLambdaSplitMassPlotsPerPtBin.add<TH2>(fmt::format("hLambdaSplitPt_vs_Cent_from_{0}_to_{1}", lambdahistvalue[i], lambdahistvalue[i + 1]).c_str(), fmt::format("Split Lambda mass vs centrality, pT from {0} to {1}", lambdahistvalue[i], lambdahistvalue[i + 1]).c_str(), {HistType::kTH2D, {{lambdaMassAxis}, {centAxis}}});
      }
      // lambdafeeddown matrices
      rFeeddownMatrices.add("hLambdaFeeddownMatrix", "hLambdaFeeddownMatrix", {HistType::kTH3D, {lambdaPtAxis, lambdaPtAxis, centAxis}});
      rFeeddownMatrices.add("hLambdaXiMinusFeeddownMatrix", "hLambdaXiMinusFeeddownMatrix", {HistType::kTH3D, {lambdaPtAxis, lambdaPtAxis, centAxis}});
      rFeeddownMatrices.add("hLambdaXiZeroFeeddownMatrix", "hLambdaXiZeroFeeddownMatrix", {HistType::kTH3D, {lambdaPtAxis, lambdaPtAxis, centAxis}});
      rFeeddownMatrices.add("hLambdaOmegaFeeddownMatrix", "hLambdaOmegaFeeddownMatrix", {HistType::kTH3D, {lambdaPtAxis, lambdaPtAxis, centAxis}});
    }
    // Adding Antilambda Histograms
    if (antiLambdaAnalysis == true) {
      // same method as in Lambda and Kzerosh above
      rPtAnalysis.add("hMassAntiLambdavsCuts", "hMassAntiLambdavsCuts", {HistType::kTH2F, {{partCutsAxis}, {k0ShortMassAxis}}});
      rPtAnalysis.add("hArmenterosPodolanskiPlotAntiLambda", "hArmenterosPodolanskiPlotAntiLambda", {HistType::kTH2F, {{armenterosasymAxis}, {armenterosQtAxis}}});
      rPtAnalysis.add("hNSigmaPosPionFromAntiLambdas", "hNSigmaPosPionFromAntiLambdas", {HistType::kTH2F, {{100, -5.f, 5.f}, {antilambdaPtAxis}}});
      rPtAnalysis.add("hNSigmaNegProtonFromAntiLambdas", "hNSigmaNegProtonFromAntiLambdas", {HistType::kTH2F, {{100, -5.f, 5.f}, {antilambdaPtAxis}}});
      rPtAnalysis.add("hAntiLambdaV0radius", "hAntiLambdaV0radius", {HistType::kTH1F, {{nBins, 0.0f, 50.0f}}});
      rPtAnalysis.add("hAntiLambdacosPA", "hAntiLambdacosPA", {HistType::kTH1F, {{nBins, 0.95f, 1.0f}}});
      rPtAnalysis.add("hAntiLambdaDCAV0Daughters", "hAntiLambdaDCAV0Daughters", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      rPtAnalysis.add("hAntiLambdaDCAPosDaughter", "hAntiLambdaDCAPosDaughter", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      rPtAnalysis.add("hAntiLambdaDCANegDaughter", "hAntiLambdaDCANegDaughter", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      for (int i = 0; i < nAntilambdaHistograms; i++) {
        pthistos::antilambdaPt[i] = rAntilambdaMassPlotsPerPtBin.add<TH2>(fmt::format("hAntiLambdaPt_vs_Cent_from_{0}_to_{1}", antilambdahistvalue[i], antilambdahistvalue[i + 1]).c_str(), fmt::format("AntiLambda mass vs centrality, pT from {0} to {1}", antilambdahistvalue[i], antilambdahistvalue[i + 1]).c_str(), {HistType::kTH2D, {{antiLambdaMassAxis}, {centAxis}}});
        pthistos::antilambdaSplit[i] = rAntilambdaSplitMassPlotsPerPtBin.add<TH2>(fmt::format("hAntiLambdaSplitPt_vs_Cent_from_{0}_to_{1}", antilambdahistvalue[i], antilambdahistvalue[i + 1]).c_str(), fmt::format("Split AntiLambda mass vs centrality, pT from {0} to {1}", antilambdahistvalue[i], antilambdahistvalue[i + 1]).c_str(), {HistType::kTH2D, {{antiLambdaMassAxis}, {centAxis}}});
      }
      // antilambdafeeddown matrices
      rFeeddownMatrices.add("hAntiLambdaFeeddownMatrix", "hAntiLambdaFeeddownMatrix", {HistType::kTH3D, {antilambdaPtAxis, antilambdaPtAxis, centAxis}});
      rFeeddownMatrices.add("hAntiLambdaXiPlusFeeddownMatrix", "hAntiLambdaXiPlusFeeddownMatrix", {HistType::kTH3D, {antilambdaPtAxis, antilambdaPtAxis, centAxis}});
      rFeeddownMatrices.add("hAntiLambdaAntiXiZeroFeeddownMatrix", "hAntiLambdaAntiXiZeroFeeddownMatrix", {HistType::kTH3D, {antilambdaPtAxis, antilambdaPtAxis, centAxis}});
      rFeeddownMatrices.add("hAntiLambdaAntiOmegaFeeddownMatrix", "hAntiLambdaAntiOmegaFeeddownMatrix", {HistType::kTH3D, {antilambdaPtAxis, antilambdaPtAxis, centAxis}});
    }

    // Particle Level Corrections
    rMCCorrections.add("hK0shBeforeEventSelectionPtSpectrum", "hK0shBeforeEventSelectionPtSpectrum", {HistType::kTH2D, {k0ShortPtAxis, centAxis}}); // not filled
    rMCCorrections.add("hLambdaBeforeEventSelectionPtSpectrum", "hLambdaBeforeEventSelectionPtSpectrum", {HistType::kTH2D, {lambdaPtAxis, centAxis}});
    rMCCorrections.add("hAntiLambdaBeforeEventSelectionPtSpectrum", "hAntiLambdaBeforeEventSelectionPtSpectrum", {HistType::kTH2D, {antilambdaPtAxis, centAxis}});
    rMCCorrections.add("hK0shAfterEventSelectionPtSpectrum", "hK0shAfterEventSelectionPtSpectrum", {HistType::kTH2D, {k0ShortPtAxis, centAxis}});
    rMCCorrections.add("hLambdaAfterEventSelectionPtSpectrum", "hLambdaAfterEventSelectionPtSpectrum", {HistType::kTH2D, {lambdaPtAxis, centAxis}});
    rMCCorrections.add("hAntiLambdaAfterEventSelectionPtSpectrum", "hAntiLambdaAfterEventSelectionPtSpectrum", {HistType::kTH2D, {antilambdaPtAxis, centAxis}});

    // Event and V0s Corrections
    rMCCorrections.add("hNEvents_Corrections", "hNEvents_Corrections", {HistType::kTH2D, {{10, 0.f, 10.f}, centAxis}});

    // Generated Level Pt Spectrums (with rapidity cut)
    rMCCorrections.add("GenParticleRapidity", "GenParticleRapidity", {HistType::kTH1F, {{nBins, -10.0f, 10.0f}}});
    rMCCorrections.add("hK0shGeneratedPtSpectrum", "hK0shGeneratedPtSpectrum", {HistType::kTH2D, {k0ShortPtAxis, centAxis}});
    rMCCorrections.add("hLambdaGeneratedPtSpectrum", "hLambdaGeneratedPtSpectrum", {HistType::kTH2D, {lambdaPtAxis, centAxis}});
    rMCCorrections.add("hAntiLambdaGeneratedPtSpectrum", "hAntiLambdaGeneratedPtSpectrum", {HistType::kTH2D, {antilambdaPtAxis, centAxis}});
    rMCCorrections.add("hXiMinusGeneratedPtSpectrum", "hXiMinusGeneratedPtSpectrum", {HistType::kTH2D, {lambdaPtAxis, centAxis}});
    rMCCorrections.add("hXiZeroGeneratedPtSpectrum", "hXiZeroGeneratedPtSpectrum", {HistType::kTH2D, {lambdaPtAxis, centAxis}});
    rMCCorrections.add("hOmegaGeneratedPtSpectrum", "hOmegaGeneratedPtSpectrum", {HistType::kTH2D, {lambdaPtAxis, centAxis}});
    rMCCorrections.add("hXiPlusGeneratedPtSpectrum", "hXiPlusGeneratedPtSpectrum", {HistType::kTH2D, {antilambdaPtAxis, centAxis}});
    rMCCorrections.add("hAntiXiZeroGeneratedPtSpectrum", "hAntiXiZeroGeneratedPtSpectrum", {HistType::kTH2D, {antilambdaPtAxis, centAxis}});
    rMCCorrections.add("hAntiOmegaGeneratedPtSpectrum", "hAntiOmegaGeneratedPtSpectrum", {HistType::kTH2D, {antilambdaPtAxis, centAxis}});
    rMCCorrections.add("hPhiGeneratedPtSpectrum", "hPhiGeneratedPtSpectrum", {HistType::kTH2D, {k0ShortPtAxis, centAxis}});
  }

  // Event selection function
  template <typename TCollision>
  bool acceptEvent(TCollision const& collision)
  {
    rPtAnalysis.fill(HIST("hNEvents"), 0.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "All");
    if (dosel8 && !collision.sel8()) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 1.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "sel 8");
    if (doNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 2.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "NoTimeFrameBorder");
    if (doNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 3.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(4, "NoITSROFrameBorder");
    if (doIsTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 4.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(5, "IsTriggerTVX");
    if (docutZVertex && !(std::abs(collision.posZ()) < cutZVertex)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 5.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(6, "cutZVertex");
    if (doisInelGt0 && !(collision.multNTracksPVeta1() > 0)) {
      // if (doisInelGt0 && !(collision.multMCNParticlesEta10() > 0)) { //CHANGE TO THIS
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 6.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(7, "isInelGt0");
    // Cut Plots
    rPtAnalysis.fill(HIST("hVertexZRec"), collision.posZ());
    return true;
  }

  // V0 selection function
  template <typename TV0, typename Track, typename TCollision>
  bool acceptV0(TV0 const& v0, Track const& posDaughterTrack, Track const& negDaughterTrack, TCollision const& collision)
  {
    rPtAnalysis.fill(HIST("hNV0s"), 0.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(1, "All V0s");
    if (doDaughterPseudorapidityCut && !(std::abs(posDaughterTrack.eta()) < etadau && std::abs(negDaughterTrack.eta()) < etadau)) { // Daughters Pseudorapidity Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNV0s"), 1.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(2, "Dau Pseudorapidity");
    if (doisNotITSAfterburner && (posDaughterTrack.isITSAfterburner() || negDaughterTrack.isITSAfterburner())) { // ITS After Burner on daughter tracks
      return false;
    }
    rPtAnalysis.fill(HIST("hNV0s"), 2.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(3, "ITS Afterburner");
    if (doitsMinHits && !(posDaughterTrack.itsNCls() >= itsMinHits && negDaughterTrack.itsNCls() >= itsMinHits)) { // Minimum hits in the ITS
      return false;
    }
    rPtAnalysis.fill(HIST("hNV0s"), 3.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(4, "ITS Min Hits");
    // Cut Plots
    rPtAnalysis.fill(HIST("hV0EtaDaughters"), v0.template posTrack_as<DaughterTracks>().eta());
    rPtAnalysis.fill(HIST("hV0EtaDaughters"), v0.template negTrack_as<DaughterTracks>().eta());
    rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlot"), v0.alpha(), v0.qtarm());
    return true;
  }

  // K0sh selection function
  template <typename TV0, typename Track, typename TCollision>
  bool acceptK0sh(TV0 const& v0, Track const& posDaughterTrack, Track const& negDaughterTrack, TCollision const& collision)
  {
    rPtAnalysis.fill(HIST("hNK0sh"), 0.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(1, "All");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 0.5, v0.mK0Short());

    if (doRapidityCut && (std::abs(v0.rapidity(0)) > rapidityCut)) { // V0 Rapidity Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 1.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(2, "Rapidity");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 1.5, v0.mK0Short());
    if (doK0shTPCPID && (std::abs(posDaughterTrack.tpcNSigmaPi()) > nSigmaTPCPion || std::abs(negDaughterTrack.tpcNSigmaPi()) > nSigmaTPCPion)) { // TPC PID for two pions
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 2.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(3, "TPC_PID");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 2.5, v0.mK0Short());
    if (doK0shcomptmasscut && ((std::abs(v0.mLambda() - o2::constants::physics::MassLambda0) < compv0masscut) || (std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) < compv0masscut))) { // Kzero competitive v0 mass cut (cut out Lambdas and Anti-Lambdas)
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 3.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(4, "Compt_Mass");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 3.5, v0.mK0Short());
    if (doK0shMaxct && (v0.v0radius() > k0shmaxct)) { // K0sh max ct
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 4.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(5, "Max_ct");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 4.5, v0.mK0Short());
    if (doK0shArmenterosCut && (v0.qtarm() < (k0shparamArmenterosCut * std::abs(v0.alpha())))) { // K0sh Armenteros Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 5.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(6, "Armenteros");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 5.5, v0.mK0Short());
    if (doK0shcosPACut && (v0.v0cosPA() < k0shSettingcosPA)) { // K0sh cosPA Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 6.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(7, "cosPA");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 6.5, v0.mK0Short());
    if (doK0shDCAdauCut && (v0.dcaV0daughters() > k0shSettingdcav0dau)) { // K0sh DCAdaughters Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 7.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(8, "DCAdau");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 7.5, v0.mK0Short());
    if (doK0shv0radiusCut && (v0.v0radius() < k0shSettingradius)) { // K0sh v0radius Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 8.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(9, "v0radius");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 8.5, v0.mK0Short());
    if (doK0shdcaposdautopv && (std::abs(v0.dcapostopv()) < k0shSettingdcapostopv)) { // K0sh DCAPosDaughterToPV Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 9.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(10, "DCAPosDautoPV");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 9.5, v0.mK0Short());
    if (doK0shdcanegdautopv && (std::abs(v0.dcanegtopv()) < k0shSettingdcanegtopv)) { // K0sh DCANegDaughterToPV Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 10.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(11, "DCANegDautoPV");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 10.5, v0.mK0Short());

    // Cut Plots
    rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlotK0sh"), v0.alpha(), v0.qtarm());
    // rPtAnalysis.fill(HIST("hNSigmaPosPionFromK0s"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tpcInnerParam());
    // rPtAnalysis.fill(HIST("hNSigmaNegPionFromK0s"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
    rPtAnalysis.fill(HIST("hK0shcosPA"), v0.v0cosPA());
    rPtAnalysis.fill(HIST("hK0shV0radius"), v0.v0radius());
    rPtAnalysis.fill(HIST("hK0shDCAV0Daughters"), v0.dcaV0daughters());
    rPtAnalysis.fill(HIST("hK0shDCAPosDaughter"), v0.dcapostopv());
    rPtAnalysis.fill(HIST("hK0shDCANegDaughter"), v0.dcanegtopv());
    rPtAnalysis.fill(HIST("V0Rapidity"), v0.rapidity(0));
    return true;
  }

  // Lambda selection function
  template <typename TV0, typename Track, typename TCollision>
  bool acceptLambda(TV0 const& v0, Track const& posDaughterTrack, Track const& negDaughterTrack, TCollision const& collision)
  {
    rPtAnalysis.fill(HIST("hNLambda"), 0.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(1, "All");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 0.5, v0.mLambda());

    if (doRapidityCut && (std::abs(v0.rapidity(1)) > rapidityCut)) { // V0 Rapidity Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 1.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(2, "Rapidity");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 1.5, v0.mLambda());
    if (doLambdaTPCPID && ((std::abs(posDaughterTrack.tpcNSigmaPr()) > nSigmaTPCProton) || (std::abs(negDaughterTrack.tpcNSigmaPi()) > nSigmaTPCPion))) { // TPC PID on daughter pion and proton for Lambda
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 2.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(3, "TPC_PID");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 2.5, v0.mLambda());
    if (doLambdacomptmasscut && ((std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < compv0masscut) || (std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) < compv0masscut))) { // Lambda competitive v0 mass cut (cut out Kaons)
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 3.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(4, "Compt_Mass");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 3.5, v0.mLambda());
    if (doLambdaMaxct && (v0.v0radius() > lambdamaxct)) { // Lambda max ct
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 4.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(5, "Max_ct");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 4.5, v0.mLambda());
    if (doLambdaArmenterosCut && v0.qtarm() < (lambdaparamArmenterosCut * std::abs(v0.alpha()))) { // Lambda Armenteros Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 5.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(6, "Armenteros");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 5.5, v0.mLambda());
    if (doLambdacosPACut && (v0.v0cosPA() < lambdaSettingcosPA)) { // Lambda cosPA Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 6.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(7, "cosPA");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 6.5, v0.mLambda());
    if (doLambdaDCAdauCut && (v0.dcaV0daughters() > lambdaSettingdcav0dau)) { // Lambda DCAdaughters Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 7.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(8, "DCAdau");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 7.5, v0.mLambda());
    if (doLambdav0radiusCut && (v0.v0radius() < lambdaSettingradius)) { // Lambda v0radius Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 8.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(9, "v0radius");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 8.5, v0.mLambda());
    if (doLambdadcaposdautopv && (std::abs(v0.dcapostopv()) < lambdaSettingdcapostopv)) { // Lambda DCAPosDaughterToPV Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 9.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(10, "DCAPosDautoPV");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 9.5, v0.mLambda());
    if (doLambdadcanegdautopv && (std::abs(v0.dcanegtopv()) < lambdaSettingdcanegtopv)) { // Lambda DCANegDaughterToPV Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 10.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(11, "DCANegDautoPV");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 10.5, v0.mLambda());

    // Cut Plots
    rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlotLambda"), v0.alpha(), v0.qtarm());
    // rPtAnalysis.fill(HIST("hNSigmaPosProtonFromLambdas"), posDaughterTrack.tpcNSigmaPr(), posDaughterTrack.tpcInnerParam());
    // rPtAnalysis.fill(HIST("hNSigmaNegPionFromLambdas"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
    rPtAnalysis.fill(HIST("hLambdacosPA"), v0.v0cosPA());
    rPtAnalysis.fill(HIST("hLambdaV0radius"), v0.v0radius());
    rPtAnalysis.fill(HIST("hLambdaDCAV0Daughters"), v0.dcaV0daughters());
    rPtAnalysis.fill(HIST("hLambdaDCAPosDaughter"), v0.dcapostopv());
    rPtAnalysis.fill(HIST("hLambdaDCANegDaughter"), v0.dcanegtopv());
    rPtAnalysis.fill(HIST("V0Rapidity"), v0.rapidity(1));
    return true;
  }

  // Antilambda selection function
  template <typename TV0, typename Track, typename TCollision>
  bool acceptAntilambda(TV0 const& v0, Track const& posDaughterTrack, Track const& negDaughterTrack, TCollision const& collision)
  {
    rPtAnalysis.fill(HIST("hNAntiLambda"), 0.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNAntiLambda"))->GetXaxis()->SetBinLabel(1, "All");
    rPtAnalysis.fill(HIST("hMassAntiLambdavsCuts"), 0.5, v0.mAntiLambda());

    if (doRapidityCut && (std::abs(v0.rapidity(2)) > rapidityCut)) { // V0 Rapidity Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntiLambda"), 1.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNAntiLambda"))->GetXaxis()->SetBinLabel(2, "Rapidity");
    rPtAnalysis.fill(HIST("hMassAntiLambdavsCuts"), 1.5, v0.mAntiLambda());
    if (doAntilambdaTPCPID && (std::abs(negDaughterTrack.tpcNSigmaPr()) > nSigmaTPCProton || std::abs(posDaughterTrack.tpcNSigmaPi()) > nSigmaTPCPion)) { // TPC PID on daughter pion and proton for AntiLambda
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntiLambda"), 2.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNAntiLambda"))->GetXaxis()->SetBinLabel(3, "TPC_PID");
    rPtAnalysis.fill(HIST("hMassAntiLambdavsCuts"), 2.5, v0.mAntiLambda());
    if (doAntilambdacomptmasscut && ((std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < compv0masscut) || (std::abs(v0.mLambda() - o2::constants::physics::MassLambda0) < compv0masscut))) { // Antilambda competitive v0 mass cut (cut out Kaons)
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntiLambda"), 3.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNAntiLambda"))->GetXaxis()->SetBinLabel(4, "Compt_Mass");
    rPtAnalysis.fill(HIST("hMassAntiLambdavsCuts"), 3.5, v0.mAntiLambda());
    if (doAntilambdaMaxct && (v0.v0radius() > antilambdamaxct)) { // Antilambda max ct
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntiLambda"), 4.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNAntiLambda"))->GetXaxis()->SetBinLabel(5, "Max_ct");
    rPtAnalysis.fill(HIST("hMassAntiLambdavsCuts"), 4.5, v0.mAntiLambda());
    if (doAntilambdaArmenterosCut && (v0.qtarm() < (antilambdaparamArmenterosCut * std::abs(v0.alpha())))) { // Antilambda Armenteros Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntiLambda"), 5.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNAntiLambda"))->GetXaxis()->SetBinLabel(6, "Armenteros");
    rPtAnalysis.fill(HIST("hMassAntiLambdavsCuts"), 5.5, v0.mAntiLambda());
    if (doAntilambdacosPACut && (v0.v0cosPA() < antilambdaSettingcosPA)) { // Antilambda cosPA Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntiLambda"), 6.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNAntiLambda"))->GetXaxis()->SetBinLabel(7, "cosPA");
    rPtAnalysis.fill(HIST("hMassAntiLambdavsCuts"), 6.5, v0.mAntiLambda());
    if (doAntilambdaDCAdauCut && (v0.dcaV0daughters() > antilambdaSettingdcav0dau)) { // Antilambda DCAdaughters Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntiLambda"), 7.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNAntiLambda"))->GetXaxis()->SetBinLabel(8, "DCAdau");
    rPtAnalysis.fill(HIST("hMassAntiLambdavsCuts"), 7.5, v0.mAntiLambda());
    if (doAntilambdav0radiusCut && (v0.v0radius() < antilambdaSettingradius)) { // Antilambda v0radius Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntiLambda"), 8.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNAntiLambda"))->GetXaxis()->SetBinLabel(9, "v0radius");
    rPtAnalysis.fill(HIST("hMassAntiLambdavsCuts"), 8.5, v0.mAntiLambda());
    if (doAntilambdadcaposdautopv && (std::abs(v0.dcapostopv()) < antilambdaSettingdcapostopv)) { // Antilambda DCAPosDaughterToPV Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntilambda"), 9.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNAntilambda"))->GetXaxis()->SetBinLabel(10, "DCAPosDautoPV");
    rPtAnalysis.fill(HIST("hMassAntilambdavsCuts"), 9.5, v0.mAntiLambda());
    if (doAntilambdadcanegdautopv && (std::abs(v0.dcanegtopv()) < antilambdaSettingdcanegtopv)) { // Antilambda DCANegDaughterToPV Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntilambda"), 10.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNAntilambda"))->GetXaxis()->SetBinLabel(11, "DCANegDautoPV");
    rPtAnalysis.fill(HIST("hMassAntilambdavsCuts"), 10.5, v0.mAntiLambda());

    // Cut plots
    rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlotAntilambda"), v0.alpha(), v0.qtarm());
    // rPtAnalysis.fill(HIST("hNSigmaPosPionFromAntilambdas"), posDaughterTrack.tpcNSigmaPr(), posDaughterTrack.tpcInnerParam());
    // rPtAnalysis.fill(HIST("hNSigmaNegProtonFromAntilambdas"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
    rPtAnalysis.fill(HIST("hAntilambdacosPA"), v0.v0cosPA());
    rPtAnalysis.fill(HIST("hAntilambdaV0radius"), v0.v0radius());
    rPtAnalysis.fill(HIST("hAntilambdaDCAV0Daughters"), v0.dcaV0daughters());
    rPtAnalysis.fill(HIST("hAntilambdaDCAPosDaughter"), v0.dcapostopv());
    rPtAnalysis.fill(HIST("hAntilambdaDCANegDaughter"), v0.dcanegtopv());
    rPtAnalysis.fill(HIST("V0Rapidity"), v0.rapidity(2));
    return true;
  }

  // V0 selection function
  template <typename TV0, typename Track, typename TCollision>
  bool acceptV0Derived(TV0 const& v0, Track const& posDaughterTrack, Track const& negDaughterTrack, TCollision const& collision)
  {
    rPtAnalysis.fill(HIST("hNV0s"), 0.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(1, "All V0s");
    if (doDaughterPseudorapidityCut && !(std::abs(v0.positiveeta()) < etadau && std::abs(v0.negativeeta()) < etadau)) { // Daughters Pseudorapidity Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNV0s"), 1.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(2, "Dau Pseudorapidity");
    if (doisNotITSAfterburner && (posDaughterTrack.isITSAfterburner() || negDaughterTrack.isITSAfterburner())) { // ITS After Burner on daughter tracks
      return false;
    }
    rPtAnalysis.fill(HIST("hNV0s"), 2.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(3, "ITS Afterburner");
    if (doitsMinHits && !(posDaughterTrack.itsNCls() >= itsMinHits && negDaughterTrack.itsNCls() >= itsMinHits)) { // Minimum hits in the ITS
      return false;
    }
    rPtAnalysis.fill(HIST("hNV0s"), 3.5, collision.centFT0M());
    rPtAnalysis.get<TH2>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(4, "ITS Min Hits");
    // Cut Plots
    rPtAnalysis.fill(HIST("hV0EtaDaughters"), v0.positiveeta());
    rPtAnalysis.fill(HIST("hV0EtaDaughters"), v0.negativeeta());
    rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlot"), v0.alpha(), v0.qtarm());
    return true;
  }

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;
  using DaughterTracksDerived = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
  o2::framework::Service<o2::framework::O2DatabasePDG> pdgDB;

  // This is the process for Generated Particles
  void genMCProcess(soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::MultMCExtras>::iterator const& mcCollision,
                    soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::PVMults, aod::McCentFT0Ms>> const& collisions,
                    aod::McParticles const& mcParticles)
  {
    // Event Efficiency, Event Split and V0 Signal Loss Corrections
    rMCCorrections.fill(HIST("hNEvents_Corrections"), 0.5, mcCollision.centFT0M()); // All Events
    if (std::abs(mcCollision.posZ()) > cutZVertex) {
      return;
    }
    if (!(mcCollision.multMCNParticlesEta10() > 0)) { // TRY TO CHANGE TO THIS
      // if (!pwglf::isINELgtNmc(mcParticles, 0, pdgDB)) {
      return;
    }
    rMCCorrections.fill(HIST("hNEvents_Corrections"), 1.5, mcCollision.centFT0M()); // Event Efficiency Denominator
    // Particles (of interest) Generated Pt Spectrum and Signal Loss Denominator Loop
    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) < rapidityCut) {
        if (mcParticle.isPhysicalPrimary()) {
          rMCCorrections.fill(HIST("GenParticleRapidity"), mcParticle.y());
          if (mcParticle.pdgCode() == kK0Short) // kzero matched
          {
            rMCCorrections.fill(HIST("hK0shGeneratedPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
          }
          if (mcParticle.pdgCode() == kLambda0) // lambda matched
          {
            rMCCorrections.fill(HIST("hLambdaGeneratedPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
          }
          if (mcParticle.pdgCode() == kLambda0Bar) // antilambda matched
          {
            rMCCorrections.fill(HIST("hAntiLambdaGeneratedPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
          }
          if (mcParticle.pdgCode() == kXiMinus) // Xi Minus matched
          {
            rMCCorrections.fill(HIST("hXiMinusGeneratedPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
          }
          if (mcParticle.pdgCode() == kXi0) // Xi Zero matched
          {
            rMCCorrections.fill(HIST("hXiZeroGeneratedPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
          }
          if (mcParticle.pdgCode() == kOmegaMinus) // Omega matched
          {
            rMCCorrections.fill(HIST("hOmegaGeneratedPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
          }
          if (mcParticle.pdgCode() == kXiPlusBar) // Xi Plus matched
          {
            rMCCorrections.fill(HIST("hXiPlusGeneratedPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
          }
          if (mcParticle.pdgCode() == -kXi0) // Anti-Xi Zero matched
          {
            rMCCorrections.fill(HIST("hAntiXiZeroGeneratedPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
          }
          if (mcParticle.pdgCode() == kOmegaPlusBar) // Anti-Omega matched
          {
            rMCCorrections.fill(HIST("hAntiOmegaGeneratedPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
          }
          if (mcParticle.pdgCode() == kPhi) // Phi
          {
            rMCCorrections.fill(HIST("hPhiGeneratedPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
          }
        }
      }
    }
    // Signal Loss Numenator Loop
    for (const auto& collision : collisions) {
      rMCCorrections.fill(HIST("hNEvents_Corrections"), 2.5, mcCollision.centFT0M()); // Number of Events Reconsctructed
      if (!acceptEvent(collision)) {                                                  // Event Selection
        return;
      }
      rMCCorrections.fill(HIST("hNEvents_Corrections"), 3.5, mcCollision.centFT0M()); // Event Split Denomimator and Event Efficiency Numenator
      for (const auto& mcParticle : mcParticles) {
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }
        if (std::abs(mcParticle.y()) > rapidityCut) {
          continue;
        }
        if (mcParticle.pdgCode() == kK0Short) // kzero matched
        {
          rMCCorrections.fill(HIST("hK0shAfterEventSelectionPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
        }
        if (mcParticle.pdgCode() == kLambda0) // lambda matched
        {
          rMCCorrections.fill(HIST("hLambdaAfterEventSelectionPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
        }
        if (mcParticle.pdgCode() == kLambda0Bar) // antilambda matched
        {
          rMCCorrections.fill(HIST("hAntiLambdaAfterEventSelectionPtSpectrum"), mcParticle.pt(), mcCollision.centFT0M());
        }
      }
    }
    // End of Signal Loss Numenator Loop
  }
  // This is the Process for the MC reconstructed Data
  void recMCProcess(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::PVMults, aod::MultsExtraMC, aod::CentFT0Ms>::iterator const& collision,
                    soa::Join<aod::McCollisions, aod::McCentFT0Ms> const& /*mcCollisions*/,
                    soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s,
                    DaughterTracks const&, // no need to define a variable for tracks, if we don't access them directly
                    aod::McParticles const& /*mcParticles*/)
  {
    // tokenise strings into individual values
    pthistos::kaonPtBins = o2::utils::Str::tokenize(kzeroSettingPtBinsString, ',');
    pthistos::lambdaPtBins = o2::utils::Str::tokenize(lambdaSettingPtBinsString, ',');
    pthistos::antilambdaPtBins = o2::utils::Str::tokenize(antilambdaSettingPtBinsString, ',');
    pthistos::kaonPtBins = o2::utils::Str::tokenize(kzeroSettingPtBinsString, ',');

    // Calculate number of histograms for each particle type
    int nKaonHistograms = pthistos::kaonPtBins.size() - 1;
    int nLambdaHistograms = pthistos::lambdaPtBins.size() - 1;
    int nAntilambdaHistograms = pthistos::antilambdaPtBins.size() - 1;

    // initialize and convert tokenized strings into vector of doubles for Pt Bin Edges
    std::vector<double> kaonptedgevalues(nKaonHistograms + 1);
    std::vector<double> lambdaptedgevalues(nLambdaHistograms + 1);
    std::vector<double> antilambdaptedgevalues(nAntilambdaHistograms + 1);

    // For centrality estimation
    const auto& mcCollision = collision.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();

    for (int i = 0; i < nKaonHistograms + 1; i++) {
      kaonptedgevalues[i] = std::stod(pthistos::kaonPtBins[i]);
    }
    for (int i = 0; i < nLambdaHistograms + 1; i++) {
      lambdaptedgevalues[i] = std::stod(pthistos::lambdaPtBins[i]);
    }
    for (int i = 0; i < nAntilambdaHistograms + 1; i++) {
      antilambdaptedgevalues[i] = std::stod(pthistos::antilambdaPtBins[i]);
    }
    if (!acceptEvent(collision)) { // Event Selection
      return;
    }
    rPtAnalysis.fill(HIST("hNRecEvents"), 0.5, mcCollision.centFT0M()); // Event Split Numenator
    for (const auto& v0 : V0s) {
      // Checking that the V0 is a true K0s/Lambdas/Antilambdas and then filling the parameter histograms and the invariant mass plots for different cuts (which are taken from namespace)
      const auto& posDaughterTrack = v0.template posTrack_as<DaughterTracks>();
      const auto& negDaughterTrack = v0.template negTrack_as<DaughterTracks>();
      if (!acceptV0(v0, posDaughterTrack, negDaughterTrack, mcCollision)) { // V0 Selections
        continue;
      }
      // kzero analysis
      if (kzeroAnalysis == true) {
        if (acceptK0sh(v0, posDaughterTrack, negDaughterTrack, mcCollision)) { // K0sh Selection
          // K0sh Signal Split Numerator Start
          for (int i = 0; i < nKaonHistograms; i++) {
            if (kaonptedgevalues[i] <= v0.pt() && v0.pt() < kaonptedgevalues[i + 1]) { // finding v0s with pt within the range of our bin edges for K0sh Splitting Numerator
              pthistos::kaonSplit[i]->Fill(v0.mK0Short(), mcCollision.centFT0M());     // filling the k0s namespace histograms for K0sh Splitting Numerator
            }
          }
          // K0sh Signla Split Numerator End
          if (v0.has_mcParticle()) {
            auto v0mcParticle = v0.mcParticle();
            if (dotruthk0sh && (v0mcParticle.pdgCode() == kK0Short)) { // kzero matched
              if (v0mcParticle.isPhysicalPrimary()) {
                for (int i = 0; i < nKaonHistograms; i++) {
                  if (kaonptedgevalues[i] <= v0.pt() && v0.pt() < kaonptedgevalues[i + 1]) { // finding v0s with pt within the range of our bin edges
                    pthistos::kaonPt[i]->Fill(v0.mK0Short(), mcCollision.centFT0M());        // filling the k0s namespace histograms
                  }
                }
              }
              if (!v0mcParticle.isPhysicalPrimary()) {
                auto v0mothers = v0mcParticle.mothers_as<aod::McParticles>(); // Get mothers
                if (!v0mothers.empty()) {
                  auto& v0mcParticleMother = v0mothers.front(); // First mother
                  rFeeddownMatrices.fill(HIST("hK0shFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt(), mcCollision.centFT0M());
                  if (v0mcParticleMother.pdgCode() == kPhi) { // Phi Mother Matched
                    rFeeddownMatrices.fill(HIST("hK0shPhiFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt(), mcCollision.centFT0M());
                  }
                }
              }
            }
          }
        }
      }
      // lambda analysis
      if (lambdaAnalysis == true) {
        if (acceptLambda(v0, posDaughterTrack, negDaughterTrack, mcCollision)) { // Lambda Selections
          // Lambda Signal Split Numerator Start
          for (int i = 0; i < nLambdaHistograms; i++) {
            if (lambdaptedgevalues[i] <= v0.pt() && v0.pt() < lambdaptedgevalues[i + 1]) {
              pthistos::lambdaSplit[i]->Fill(v0.mLambda(), mcCollision.centFT0M());
            }
          }
          // Lambda Signal Split Numerator End
          if (v0.has_mcParticle()) {
            auto v0mcParticle = v0.mcParticle();
            if (dotruthLambda && (v0mcParticle.pdgCode() == kLambda0)) { // lambda matched
              if (v0mcParticle.isPhysicalPrimary()) {
                for (int i = 0; i < nLambdaHistograms; i++) {
                  if (lambdaptedgevalues[i] <= v0.pt() && v0.pt() < lambdaptedgevalues[i + 1]) {
                    pthistos::lambdaPt[i]->Fill(v0.mLambda(), mcCollision.centFT0M());
                  }
                }
              }
              if (!v0mcParticle.isPhysicalPrimary()) {
                auto v0mothers = v0mcParticle.mothers_as<aod::McParticles>(); // Get mothers
                if (!v0mothers.empty()) {
                  auto& v0mcParticleMother = v0mothers.front(); // First mother
                  rFeeddownMatrices.fill(HIST("hLambdaFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt(), mcCollision.centFT0M());
                  if (v0mcParticleMother.pdgCode() == kXiMinus) { // Xi Minus Mother Matched
                    rFeeddownMatrices.fill(HIST("hLambdaXiMinusFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt(), mcCollision.centFT0M());
                  }
                  if (v0mcParticleMother.pdgCode() == kXi0) { // Xi Zero Mother Matched
                    rFeeddownMatrices.fill(HIST("hLambdaXiZeroFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt(), mcCollision.centFT0M());
                  }
                  if (v0mcParticleMother.pdgCode() == kOmegaMinus) { // Omega Mother Matched
                    rFeeddownMatrices.fill(HIST("hLambdaOmegaFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt(), mcCollision.centFT0M());
                  }
                }
              }
            }
          }
        }
      }
      // antilambda analysis
      if (antiLambdaAnalysis == true) {
        if (acceptAntilambda(v0, posDaughterTrack, negDaughterTrack, mcCollision)) { // Antilambda Selections
          // Antilambda Signal Split Numerator End
          for (int i = 0; i < nAntilambdaHistograms; i++) {
            if (antilambdaptedgevalues[i] <= v0.pt() && v0.pt() < antilambdaptedgevalues[i + 1]) {
              pthistos::antilambdaSplit[i]->Fill(v0.mAntiLambda(), mcCollision.centFT0M());
            }
          }
          // Antilambda Signal Split Numerator End
          if (v0.has_mcParticle()) {
            auto v0mcParticle = v0.mcParticle();
            if (dotruthAntilambda && (v0mcParticle.pdgCode() == kLambda0Bar)) { // antilambda matched
              if (v0mcParticle.isPhysicalPrimary()) {
                for (int i = 0; i < nAntilambdaHistograms; i++) {
                  if (antilambdaptedgevalues[i] <= v0.pt() && v0.pt() < antilambdaptedgevalues[i + 1]) {
                    pthistos::antilambdaPt[i]->Fill(v0.mAntiLambda(), mcCollision.centFT0M());
                  }
                }
              }
              if (!v0mcParticle.isPhysicalPrimary()) {
                auto v0mothers = v0mcParticle.mothers_as<aod::McParticles>(); // Get mothers
                if (!v0mothers.empty()) {
                  auto& v0mcParticleMother = v0mothers.front(); // First mother
                  rFeeddownMatrices.fill(HIST("hAntiLambdaFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt(), mcCollision.centFT0M());
                  if (v0mcParticleMother.pdgCode() == kXiPlusBar) { // Xi Plus Mother Matched
                    rFeeddownMatrices.fill(HIST("hAntiLambdaXiPlusFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt(), mcCollision.centFT0M());
                  }
                  if (v0mcParticleMother.pdgCode() == -kXi0) { // Anti-Xi Zero Mother Matched
                    rFeeddownMatrices.fill(HIST("hAntiLambdaAntiXiZeroFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt(), mcCollision.centFT0M());
                  }
                  if (v0mcParticleMother.pdgCode() == kOmegaPlusBar) { // Anti-Omega (minus) Mother Matched
                    rFeeddownMatrices.fill(HIST("hAntiLambdaAntiOmegaFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt(), mcCollision.centFT0M());
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
  void dataProcess(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::CentFT0Ms, aod::MultsExtra /*,aod::CentNGlobals*/>::iterator const& collision,
                   aod::V0Datas const& V0s,
                   DaughterTracks const&)
  {
    // tokenise strings into individual values
    pthistos::kaonPtBins = o2::utils::Str::tokenize(kzeroSettingPtBinsString, ',');
    pthistos::lambdaPtBins = o2::utils::Str::tokenize(lambdaSettingPtBinsString, ',');
    pthistos::antilambdaPtBins = o2::utils::Str::tokenize(antilambdaSettingPtBinsString, ',');
    pthistos::kaonPtBins = o2::utils::Str::tokenize(kzeroSettingPtBinsString, ',');

    // Calculate number of histograms for each particle type
    int nKaonHistograms = pthistos::kaonPtBins.size() - 1;
    int nLambdaHistograms = pthistos::lambdaPtBins.size() - 1;
    int nAntilambdaHistograms = pthistos::antilambdaPtBins.size() - 1;

    // initialize and convert tokenized strings into vector of doubles for Pt Bin Edges
    std::vector<double> kaonptedgevalues(nKaonHistograms + 1);
    std::vector<double> lambdaptedgevalues(nLambdaHistograms + 1);
    std::vector<double> antilambdaptedgevalues(nAntilambdaHistograms + 1);

    for (int i = 0; i < nKaonHistograms + 1; i++) {
      kaonptedgevalues[i] = std::stod(pthistos::kaonPtBins[i]);
    }
    for (int i = 0; i < nLambdaHistograms + 1; i++) {
      lambdaptedgevalues[i] = std::stod(pthistos::lambdaPtBins[i]);
    }
    for (int i = 0; i < nAntilambdaHistograms + 1; i++) {
      antilambdaptedgevalues[i] = std::stod(pthistos::antilambdaPtBins[i]);
    }

    // if (!acceptEvent(collision)) { // Event Selection
    //   return;
    // }
    rPtAnalysis.fill(HIST("hNRecEvents"), 0.5, collision.centFT0M()); // Number of recorded events
    for (const auto& v0 : V0s) {
      // Checking that the V0 is a true K0s/Lambdas/Antilambdas and then filling the parameter histograms and the invariant mass plots for different cuts (which are taken from namespace)
      const auto& posDaughterTrack = v0.template posTrack_as<DaughterTracks>();
      const auto& negDaughterTrack = v0.template negTrack_as<DaughterTracks>();
      if (!acceptV0(v0, posDaughterTrack, negDaughterTrack, collision)) { // V0 Selection
        continue;
      }
      // kzero analysis
      if (kzeroAnalysis == true) {
        if (acceptK0sh(v0, posDaughterTrack, negDaughterTrack, collision)) { // K0sh Selection
          for (int i = 0; i < nKaonHistograms; i++) {
            if (kaonptedgevalues[i] <= v0.pt() && v0.pt() < kaonptedgevalues[i + 1]) { // finding v0s with pt within the range of our bin edges
              pthistos::kaonPt[i]->Fill(v0.mK0Short(), collision.centFT0M());          // filling the k0s namespace histograms
            }
          }
        }
      }
      // lambda analysis
      if (lambdaAnalysis == true) {
        if (acceptLambda(v0, posDaughterTrack, negDaughterTrack, collision)) { // Lambda Selection
          for (int i = 0; i < nLambdaHistograms; i++) {
            if (lambdaptedgevalues[i] <= v0.pt() && v0.pt() < lambdaptedgevalues[i + 1]) {
              pthistos::lambdaPt[i]->Fill(v0.mLambda(), collision.centFT0M());
            }
          }
        }
      }
      // anti-lambda analysis
      if (antiLambdaAnalysis == true) {
        if (acceptAntilambda(v0, posDaughterTrack, negDaughterTrack, collision)) { // Antilambda Selection
          for (int i = 0; i < nAntilambdaHistograms; i++) {
            if (lambdaptedgevalues[i] <= v0.pt() && v0.pt() < lambdaptedgevalues[i + 1]) {
              pthistos::antilambdaPt[i]->Fill(v0.mAntiLambda(), collision.centFT0M());
            }
          }
        }
      }
    }
  }
  // This is the process for Real Derived Data
  void dataProcessDerived(soa::Join<aod::StraCollisions, aod::StraEvSels, aod::StraCents>::iterator const& collision,
                          soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras> const& V0s,
                          DaughterTracksDerived const&)
  {
    // tokenise strings into individual values
    pthistos::kaonPtBins = o2::utils::Str::tokenize(kzeroSettingPtBinsString, ',');
    pthistos::lambdaPtBins = o2::utils::Str::tokenize(lambdaSettingPtBinsString, ',');
    pthistos::antilambdaPtBins = o2::utils::Str::tokenize(antilambdaSettingPtBinsString, ',');
    pthistos::kaonPtBins = o2::utils::Str::tokenize(kzeroSettingPtBinsString, ',');

    // Calculate number of histograms for each particle type
    int nKaonHistograms = pthistos::kaonPtBins.size() - 1;
    int nLambdaHistograms = pthistos::lambdaPtBins.size() - 1;
    int nAntilambdaHistograms = pthistos::antilambdaPtBins.size() - 1;

    // initialize and convert tokenized strings into vector of doubles for Pt Bin Edges
    std::vector<double> kaonptedgevalues(nKaonHistograms + 1);
    std::vector<double> lambdaptedgevalues(nLambdaHistograms + 1);
    std::vector<double> antilambdaptedgevalues(nAntilambdaHistograms + 1);

    for (int i = 0; i < nKaonHistograms + 1; i++) {
      kaonptedgevalues[i] = std::stod(pthistos::kaonPtBins[i]);
    }
    for (int i = 0; i < nLambdaHistograms + 1; i++) {
      lambdaptedgevalues[i] = std::stod(pthistos::lambdaPtBins[i]);
    }
    for (int i = 0; i < nAntilambdaHistograms + 1; i++) {
      antilambdaptedgevalues[i] = std::stod(pthistos::antilambdaPtBins[i]);
    }

    // if (!acceptEvent(collision)) { // Event Selection
    //   return;
    // }
    rPtAnalysis.fill(HIST("hNRecEvents"), 0.5, collision.centFT0M()); // Number of recorded events
    for (const auto& v0 : V0s) {
      // Checking that the V0 is a true K0s/Lambdas/Antilambdas and then filling the parameter histograms and the invariant mass plots for different cuts (which are taken from namespace)
      const auto& posDaughterTrack = v0.template posTrackExtra_as<DaughterTracksDerived>(); // Positive Daughter track
      const auto& negDaughterTrack = v0.template negTrackExtra_as<DaughterTracksDerived>(); // Negative Daughter track
      if (!acceptV0Derived(v0, posDaughterTrack, negDaughterTrack, collision)) {            // V0 Selection
        continue;
      }
      // kzero analysis
      if (kzeroAnalysis == true) {
        if (acceptK0sh(v0, posDaughterTrack, negDaughterTrack, collision)) { // K0sh Selection
          for (int i = 0; i < nKaonHistograms; i++) {
            if (kaonptedgevalues[i] <= v0.pt() && v0.pt() < kaonptedgevalues[i + 1]) { // finding v0s with pt within the range of our bin edges
              pthistos::kaonPt[i]->Fill(v0.mK0Short(), collision.centFT0M());          // filling the k0s namespace histograms
            }
          }
        }
      }
      // lambda analysis
      if (lambdaAnalysis == true) {
        if (acceptLambda(v0, posDaughterTrack, negDaughterTrack, collision)) { // Lambda Selection
          for (int i = 0; i < nLambdaHistograms; i++) {
            if (lambdaptedgevalues[i] <= v0.pt() && v0.pt() < lambdaptedgevalues[i + 1]) {
              pthistos::lambdaPt[i]->Fill(v0.mLambda(), collision.centFT0M());
            }
          }
        }
      }
      // anti-lambda analysis
      if (antiLambdaAnalysis == true) {
        if (acceptAntilambda(v0, posDaughterTrack, negDaughterTrack, collision)) { // Antilambda Selection
          for (int i = 0; i < nAntilambdaHistograms; i++) {
            if (antilambdaptedgevalues[i] <= v0.pt() && v0.pt() < antilambdaptedgevalues[i + 1]) {
              pthistos::antilambdaPt[i]->Fill(v0.mAntiLambda(), collision.centFT0M());
            }
          }
        }
      }
    }
  }
  void recMCProcessDerived(soa::Join<aod::StraCollisions, aod::StraEvSels, aod::StraCollLabels, aod::StraCents>::iterator const& collision,
                           // To add McCentFT0Ms
                           soa::Join<aod::V0CollRefs, aod::V0MCCores, aod::V0Cores, aod::V0Extras, aod::V0CoreMCLabels, aod::V0MCMothers> const& V0s,
                           DaughterTracksDerived const&)
  {
    // tokenise strings into individual values
    pthistos::kaonPtBins = o2::utils::Str::tokenize(kzeroSettingPtBinsString, ',');
    pthistos::lambdaPtBins = o2::utils::Str::tokenize(lambdaSettingPtBinsString, ',');
    pthistos::antilambdaPtBins = o2::utils::Str::tokenize(antilambdaSettingPtBinsString, ',');
    pthistos::kaonPtBins = o2::utils::Str::tokenize(kzeroSettingPtBinsString, ',');

    // Calculate number of histograms for each particle type
    int nKaonHistograms = pthistos::kaonPtBins.size() - 1;
    int nLambdaHistograms = pthistos::lambdaPtBins.size() - 1;
    int nAntilambdaHistograms = pthistos::antilambdaPtBins.size() - 1;

    // initialize and convert tokenized strings into vector of doubles for Pt Bin Edges
    std::vector<double> kaonptedgevalues(nKaonHistograms + 1);
    std::vector<double> lambdaptedgevalues(nLambdaHistograms + 1);
    std::vector<double> antilambdaptedgevalues(nAntilambdaHistograms + 1);

    for (int i = 0; i < nKaonHistograms + 1; i++) {
      kaonptedgevalues[i] = std::stod(pthistos::kaonPtBins[i]);
    }
    for (int i = 0; i < nLambdaHistograms + 1; i++) {
      lambdaptedgevalues[i] = std::stod(pthistos::lambdaPtBins[i]);
    }
    for (int i = 0; i < nAntilambdaHistograms + 1; i++) {
      antilambdaptedgevalues[i] = std::stod(pthistos::antilambdaPtBins[i]);
    }
    // if (!acceptEvent(collision)) { // Event Selection
    //   return;
    // }
    rPtAnalysis.fill(HIST("hNRecEvents"), 0.5, collision.centFT0M()); // Event Split Numenator
    for (const auto& v0 : V0s) {
      // Checking that the V0 is a true K0s/Lambdas/Antilambdas and then filling the parameter histograms and the invariant mass plots for different cuts (which are taken from namespace)
      const auto& posDaughterTrack = v0.template posTrackExtra_as<DaughterTracksDerived>(); // Positive Daughter track
      const auto& negDaughterTrack = v0.template negTrackExtra_as<DaughterTracksDerived>(); // Negative Daughter track
      if (!acceptV0Derived(v0, posDaughterTrack, negDaughterTrack, collision)) {            // V0 Selections
        continue;
      }
      // kzero analysis
      if (kzeroAnalysis == true) {
        if (acceptK0sh(v0, posDaughterTrack, negDaughterTrack, collision)) { // K0sh Selection
          // K0sh Signal Split Numerator Start
          for (int i = 0; i < nKaonHistograms; i++) {
            if (kaonptedgevalues[i] <= v0.ptMC() && v0.ptMC() < kaonptedgevalues[i + 1]) { // finding v0s with pt within the range of our bin edges for K0sh Splitting Numerator
              pthistos::kaonSplit[i]->Fill(v0.mK0Short(), collision.centFT0M());           // filling the k0s namespace histograms for K0sh Splitting Numerator
            }
          }
          // K0sh SignaL Split Numerator End
          if (v0.has_v0MCCore()) {
            auto v0mcParticle = v0.v0MCCore_as<aod::V0MCCores>();
            if (dotruthk0sh && (v0mcParticle.pdgCode() == kK0Short)) { // kzero matched
              if (v0mcParticle.isPhysicalPrimary()) {
                for (int i = 0; i < nKaonHistograms; i++) {
                  if (kaonptedgevalues[i] <= v0.ptMC() && v0.ptMC() < kaonptedgevalues[i + 1]) { // finding v0s with pt within the range of our bin edges
                    pthistos::kaonPt[i]->Fill(v0.mK0Short(), collision.centFT0M());              // filling the k0s namespace histograms
                  }
                }
              }
              if (!v0mcParticle.isPhysicalPrimary()) {
                auto v0mother = v0.motherMCPart(); // Get mothers
                rFeeddownMatrices.fill(HIST("hK0shFeeddownMatrix"), v0mcParticle.ptMC(), std::hypot(v0mother.px(), v0mother.py()), collision.centFT0M());
                if (v0mother.pdgCode() == kPhi) { // Phi Mother Matched
                  rFeeddownMatrices.fill(HIST("hK0shFeeddownMatrix"), v0mcParticle.ptMC(), std::hypot(v0mother.px(), v0mother.py()), collision.centFT0M());
                }
              }
            }
          }
        }
      }
      // lambda analysis
      if (lambdaAnalysis == true) {
        if (acceptLambda(v0, posDaughterTrack, negDaughterTrack, collision)) { // Lambda Selections
          // Lambda Signal Split Numerator Start
          for (int i = 0; i < nLambdaHistograms; i++) {
            if (lambdaptedgevalues[i] <= v0.ptMC() && v0.ptMC() < lambdaptedgevalues[i + 1]) {
              pthistos::lambdaSplit[i]->Fill(v0.mLambda(), collision.centFT0M());
            }
          }
          // Lambda Signal Split Numerator End
          if (v0.has_v0MCCore()) {
            auto v0mcParticle = v0.v0MCCore_as<aod::V0MCCores>();
            if (dotruthLambda && (v0mcParticle.pdgCode() == kLambda0)) { // lambda matched
              if (v0mcParticle.isPhysicalPrimary()) {
                for (int i = 0; i < nLambdaHistograms; i++) {
                  if (lambdaptedgevalues[i] <= v0.ptMC() && v0.ptMC() < lambdaptedgevalues[i + 1]) {
                    pthistos::lambdaPt[i]->Fill(v0.mLambda(), collision.centFT0M());
                  }
                }
              }
              if (!v0mcParticle.isPhysicalPrimary()) {
                auto v0mother = v0.motherMCPart(); // Get mothers
                rFeeddownMatrices.fill(HIST("hLambdaFeeddownMatrix"), v0mcParticle.ptMC(), std::hypot(v0mother.px(), v0mother.py()), collision.centFT0M());
                if (v0mother.pdgCode() == kXiMinus) { // Xi Minus Mother Matched
                  rFeeddownMatrices.fill(HIST("hLambdaXiMinusFeeddownMatrix"), v0mcParticle.ptMC(), std::hypot(v0mother.px(), v0mother.py()), collision.centFT0M());
                }
                if (v0mother.pdgCode() == kXi0) { // Xi Zero Mother Matched
                  rFeeddownMatrices.fill(HIST("hLambdaXiZeroFeeddownMatrix"), v0mcParticle.ptMC(), std::hypot(v0mother.px(), v0mother.py()), collision.centFT0M());
                }
                if (v0mother.pdgCode() == kOmegaMinus) { // Omega Mother Matched
                  rFeeddownMatrices.fill(HIST("hLambdaOmegaFeeddownMatrix"), v0mcParticle.ptMC(), std::hypot(v0mother.px(), v0mother.py()), collision.centFT0M());
                }
              }
            }
          }
        }
      }
      // antilambda analysis
      if (antiLambdaAnalysis == true) {
        if (acceptAntilambda(v0, posDaughterTrack, negDaughterTrack, collision)) { // Antilambda Selections
          // Antilambda Signal Split Numerator End
          for (int i = 0; i < nAntilambdaHistograms; i++) {
            if (antilambdaptedgevalues[i] <= v0.ptMC() && v0.ptMC() < antilambdaptedgevalues[i + 1]) {
              pthistos::antilambdaSplit[i]->Fill(v0.mAntiLambda(), collision.centFT0M());
            }
          }
          // Antilambda Signal Split Numerator End
          if (v0.has_v0MCCore()) {
            auto v0mcParticle = v0.v0MCCore_as<aod::V0MCCores>();
            if (dotruthAntilambda && (v0mcParticle.pdgCode() == kLambda0Bar)) { // antilambda matched
              if (v0mcParticle.isPhysicalPrimary()) {
                for (int i = 0; i < nAntilambdaHistograms; i++) {
                  if (antilambdaptedgevalues[i] <= v0.ptMC() && v0.ptMC() < antilambdaptedgevalues[i + 1]) {
                    pthistos::antilambdaPt[i]->Fill(v0.mAntiLambda(), collision.centFT0M());
                  }
                }
              }
              if (!v0mcParticle.isPhysicalPrimary()) {
                auto v0mother = v0.motherMCPart(); // Get mothers
                rFeeddownMatrices.fill(HIST("hAntiLambdaFeeddownMatrix"), v0mcParticle.ptMC(), std::hypot(v0mother.px(), v0mother.py()), collision.centFT0M());
                if (v0mother.pdgCode() == kXiPlusBar) { // Xi Plus Mother Matched
                  rFeeddownMatrices.fill(HIST("hAntiLambdaXiPlusFeeddownMatrix"), v0mcParticle.ptMC(), std::hypot(v0mother.px(), v0mother.py()), collision.centFT0M());
                }
                if (v0mother.pdgCode() == -kXi0) { // Anti-Xi Zero Mother Matched
                  rFeeddownMatrices.fill(HIST("hAntiLambdaAntiXiZeroFeeddownMatrix"), v0mcParticle.ptMC(), std::hypot(v0mother.px(), v0mother.py()), collision.centFT0M());
                }
                if (v0mother.pdgCode() == kOmegaPlusBar) { // Anti-Omega (minus) Mother Matched
                  rFeeddownMatrices.fill(HIST("hAntiLambdaAntiOmegaFeeddownMatrix"), v0mcParticle.ptMC(), std::hypot(v0mother.px(), v0mother.py()), collision.centFT0M());
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
  PROCESS_SWITCH(V0PtInvMassPlots, dataProcess, "Process Run 3 Data,", false);
  // PROCESS_SWITCH(V0PtInvMassPlots, genMCProcessDerived, "Process Run 3 MC Generated", false);
  PROCESS_SWITCH(V0PtInvMassPlots, recMCProcessDerived, "Process Run 3 MC Reconstructed", false);
  PROCESS_SWITCH(V0PtInvMassPlots, dataProcessDerived, "Process Run 3 Data,", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<V0PtInvMassPlots>(cfgc)};
}
