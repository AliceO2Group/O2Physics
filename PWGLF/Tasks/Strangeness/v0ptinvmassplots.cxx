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

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

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
using namespace o2::constants::physics;

struct V0PtInvMassPlots {
  // Histogram Registries
  HistogramRegistry rPtAnalysis{"PtAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rKaonshMassPlotsPerPtBin{"KaonshMassPlotsPerPtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rLambdaMassPlotsPerPtBin{"LambdaMassPlotsPerPtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rAntilambdaMassPlotsPerPtBin{"AntilambdaMassPlotsPerPtBin", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rFeeddownMatrices{"FeeddownMatrices", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry rMCCorrections{"MCCorrections", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

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
  Configurable<bool> doisITSAfterburner{"doisITSAfterburner", true, "Enable ITS Afterburner"};
  Configurable<bool> doitsMinHits{"doitsMinHits", true, "Enable ITS Minimum hits"};

  // Configurables switches for K0sh selection
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
  Configurable<bool> doAntilambdaTPCPID{"doAntilambdaTPCPID", true, "Enable Antilambda TPC PID"};
  Configurable<bool> doAntilambdacomptmasscut{"doAntilambdacomptmasscut", true, "Enable Antilambda Competitive V0 Mass Cut"};
  Configurable<bool> doAntilambdaMaxct{"doAntilambdaMaxct", true, "Enable Antilambda Max ct Cut"};
  Configurable<bool> doAntilambdaArmenterosCut{"doAntilambdaArmenterosCut", true, "Enable Antilambda Armenteros Cut"};
  Configurable<bool> doAntilambdacosPACut{"doAntilambdacosPACut", true, "Enable Antilambda cosPA Topological Cut"};
  Configurable<bool> doAntilambdaDCAdauCut{"doAntilambdaDCAdauCut", true, "Enable Antilambda DCA daughters Topological Cut"};
  Configurable<bool> doAntilambdav0radiusCut{"doAntilambdav0radiusCut", true, "Enable Antilambda v0radius Topological Cut"};
  Configurable<bool> doAntilambdadcaposdautopv{"doAntilambdadcaposdautopv", true, "Enable Antilambda DCA pos daughter to PV Topological Cut"};
  Configurable<bool> doAntilambdadcanegdautopv{"doAntilambdadcanegdautopv", true, "Enable Antilambda DCA neg daughter to PV Topological Cut"};

  // Configurable Kaonsh Cuts (best cuts determined by v0topologicalcuts task)
  Configurable<float> kaonshSettingdcav0dau{"kaonshSettingdcav0dau", 0.3, "DCA V0 Daughters"};
  Configurable<float> kaonshSettingdcapostopv{"kaonshSettingdcapostopv", 0.05, "DCA Pos To PV"};
  Configurable<float> kaonshSettingdcanegtopv{"kaonshSettingdcanegtopv", 0.05, "DCA Neg To PV"};
  Configurable<double> kaonshSettingcosPA{"kaonshSettingcosPA", 0.98, "V0 CosPA"}; // double -> N.B. dcos(x)/dx = 0 at x=0
  Configurable<float> kaonshSettingradius{"kaonshSettingradius", 0.50, "v0radius"};
  Configurable<float> kaonshmaxct{"kaonshmaxct", 20.00, "K0sh maximum ct value"};
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
    AxisSpec armenterosasymAxis = {nBinsArmenteros, -1.f, 1.f, "#it{p}^{+}_{||}-#it{p}^{-}_{||}/#it{p}^{+}_{||}+#it{p}^{-}_{||}"};
    AxisSpec vertexZAxis = {nBins, -11.0f, 11.0f, "vrtx_{Z} [cm]"};
    AxisSpec partCutsAxis{10, 0.0f, 10.0f, "Cut index"};

    std::vector<std::string> kaonhistvalue(nmaxHistograms + 1);
    std::vector<std::string> lambdahistvalue(nmaxHistograms + 1);
    std::vector<std::string> antilambdahistvalue(nmaxHistograms + 1);
    // K0short Histogram Pt Bin Edges
    for (int i = 0; i < nmaxHistograms + 1; i++) {     // Histos won't accept "." character so converting it to "_"
      std::string kaonptbin = pthistos::kaonPtBins[i]; // getting the value of the bin edge
      size_t pos = kaonptbin.find(".");                // finding the "." character
      kaonptbin[pos] = '_';                            // changing the "." character of the string-value to a "_"
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

    // General Plots
    rPtAnalysis.add("hNEvents", "hNEvents", {HistType::kTH1D, {{10, 0.f, 10.f}}});
    rPtAnalysis.add("hNRecEvents_Data", "hNRecEvents_Data", {HistType::kTH1D, {{1, 0.f, 1.f}}});
    rPtAnalysis.add("hNV0s", "hNV0s", {HistType::kTH1D, {{10, 0.f, 10.f}}});
    rPtAnalysis.add("hNK0sh", "hNK0sh", {HistType::kTH1D, {{10, 0.f, 10.f}}});
    rPtAnalysis.add("hNLambda", "hNLambda", {HistType::kTH1D, {{10, 0.f, 10.f}}});
    rPtAnalysis.add("hNAntilambda", "hNAntilambda", {HistType::kTH1D, {{10, 0.f, 10.f}}});

    rPtAnalysis.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rPtAnalysis.add("hArmenterosPodolanskiPlot", "hArmenterosPodolanskiPlot", {HistType::kTH2F, {{armenterosasymAxis}, {armenterosQtAxis}}});
    rPtAnalysis.add("hV0EtaDaughters", "hV0EtaDaughters", {HistType::kTH1F, {{nBins, -1.2f, 1.2f}}});
    rPtAnalysis.add("V0Rapidity", "V0Rapidity", {HistType::kTH1F, {{nBins, -10.0f, 10.0f}}});

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
      for (int i = 0; i < nmaxHistograms; i++) {
        pthistos::kaonPt[i] = rKaonshMassPlotsPerPtBin.add<TH1>(fmt::format("hPt_from_{0}_to_{1}", kaonhistvalue[i], kaonhistvalue[i + 1]).c_str(), fmt::format("hPt_from_{0}_to_{1}", kaonhistvalue[i], kaonhistvalue[i + 1]).c_str(), {HistType::kTH1D, {{k0ShortMassAxis}}});
      }
      rFeeddownMatrices.add("hK0shFeeddownMatrix", "hK0shFeeddownMatrix", {HistType::kTH2F, {{k0ShortPtAxis}, {k0ShortPtAxis}}});
      rFeeddownMatrices.add("hK0shPhiFeeddownMatrix", "hK0shPhiFeeddownMatrix", {HistType::kTH2F, {{k0ShortPtAxis}, {k0ShortPtAxis}}});
    }
    // Adding Lambda Histograms
    if (lambdaAnalysis == true) {
      // same method as in Kzerosh above
      rPtAnalysis.add("hMassLambdavsCuts", "hMassLambdavsCuts", {HistType::kTH2F, {{partCutsAxis}, {k0ShortMassAxis}}});
      rPtAnalysis.add("hArmenterosPodolanskiPlotLambda", "hArmenterosPodolanskiPlotLambda", {HistType::kTH2F, {{armenterosasymAxis}, {armenterosQtAxis}}});
      rPtAnalysis.add("hLambdaAlphaTestPtSpectrum", "hLambdaAlphaTestPtSpectrum", {HistType::kTH1F, {lambdaPtAxis}});
      rPtAnalysis.add("hNSigmaPosProtonFromLambdas", "hNSigmaPosProtonFromLambdas", {HistType::kTH2F, {{100, -5.f, 5.f}, {lambdaPtAxis}}});
      rPtAnalysis.add("hNSigmaNegPionFromLambdas", "hNSigmaNegPionFromLambdas", {HistType::kTH2F, {{100, -5.f, 5.f}, {lambdaPtAxis}}});
      rPtAnalysis.add("hLambdaV0radius", "hLambdaV0radius", {HistType::kTH1F, {{nBins, 0.0f, 50.0f}}});
      rPtAnalysis.add("hLambdacosPA", "hLambdacosPA", {HistType::kTH1F, {{nBins, 0.95f, 1.0f}}});
      rPtAnalysis.add("hLambdaDCAV0Daughters", "hLambdaDCAV0Daughters", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      rPtAnalysis.add("hLambdaDCAPosDaughter", "hLambdaDCAPosDaughter", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      rPtAnalysis.add("hLambdaDCANegDaughter", "hLambdaDCANegDaughter", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      for (int i = 0; i < nmaxHistograms; i++) {
        pthistos::lambdaPt[i] = rLambdaMassPlotsPerPtBin.add<TH1>(fmt::format("hPt_from_{0}_to_{1}", lambdahistvalue[i], lambdahistvalue[i + 1]).c_str(), fmt::format("hPt_from_{0}_to_{1}", lambdahistvalue[i], lambdahistvalue[i + 1]).c_str(), {HistType::kTH1D, {{lambdaMassAxis}}});
      }
      // lambdafeeddown matrices
      rFeeddownMatrices.add("hLambdaFeeddownMatrix", "hLambdaFeeddownMatrix", {HistType::kTH2F, {{lambdaPtAxis}, {lambdaPtAxis}}});
      rFeeddownMatrices.add("hLambdaXiMinusFeeddownMatrix", "hLambdaXiMinusFeeddownMatrix", {HistType::kTH2F, {{lambdaPtAxis}, {lambdaPtAxis}}});
      rFeeddownMatrices.add("hLambdaXiZeroFeeddownMatrix", "hLambdaXiZeroFeeddownMatrix", {HistType::kTH2F, {{lambdaPtAxis}, {lambdaPtAxis}}});
      rFeeddownMatrices.add("hLambdaOmegaFeeddownMatrix", "hLambdaOmegaFeeddownMatrix", {HistType::kTH2F, {{lambdaPtAxis}, {lambdaPtAxis}}});
    }
    // Adding Antilambda Histograms
    if (antiLambdaAnalysis == true) {
      // same method as in Lambda and Kzerosh above
      rPtAnalysis.add("hMassAntilambdavsCuts", "hMassAntilambdavsCuts", {HistType::kTH2F, {{partCutsAxis}, {k0ShortMassAxis}}});
      rPtAnalysis.add("hArmenterosPodolanskiPlotAntilambda", "hArmenterosPodolanskiPlotAntilambda", {HistType::kTH2F, {{armenterosasymAxis}, {armenterosQtAxis}}});
      rPtAnalysis.add("hAntilambdaAlphaTestPtSpectrum", "hAntilambdaAlphaTestPtSpectrum", {HistType::kTH1F, {antilambdaPtAxis}});
      rPtAnalysis.add("hNSigmaPosPionFromAntilambdas", "hNSigmaPosPionFromAntilambdas", {HistType::kTH2F, {{100, -5.f, 5.f}, {antilambdaPtAxis}}});
      rPtAnalysis.add("hNSigmaNegProtonFromAntilambdas", "hNSigmaNegProtonFromAntilambdas", {HistType::kTH2F, {{100, -5.f, 5.f}, {antilambdaPtAxis}}});
      rPtAnalysis.add("hAntilambdaV0radius", "hAntilambdaV0radius", {HistType::kTH1F, {{nBins, 0.0f, 50.0f}}});
      rPtAnalysis.add("hAntilambdacosPA", "hAntilambdacosPA", {HistType::kTH1F, {{nBins, 0.95f, 1.0f}}});
      rPtAnalysis.add("hAntilambdaDCAV0Daughters", "hAntilambdaDCAV0Daughters", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      rPtAnalysis.add("hAntilambdaDCAPosDaughter", "hAntilambdaDCAPosDaughter", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      rPtAnalysis.add("hAntilambdaDCANegDaughter", "hAntilambdaDCANegDaughter", {HistType::kTH1F, {{nBins, 0.0f, 2.2f}}});
      for (int i = 0; i < nmaxHistograms; i++) {
        pthistos::antilambdaPt[i] = rAntilambdaMassPlotsPerPtBin.add<TH1>(fmt::format("hPt_from_{0}_to_{1}", antilambdahistvalue[i], antilambdahistvalue[i + 1]).c_str(), fmt::format("hPt_from_{0}_to_{1}", antilambdahistvalue[i], antilambdahistvalue[i + 1]).c_str(), {HistType::kTH1D, {{antiLambdaMassAxis}}});
      }
      // antilambdafeeddown matrices
      rFeeddownMatrices.add("hAntiLambdaFeeddownMatrix", "hAntiLambdaFeeddownMatrix", {HistType::kTH2F, {{antilambdaPtAxis}, {antilambdaPtAxis}}});
      rFeeddownMatrices.add("hAntiLambdaXiPlusFeeddownMatrix", "hAntiLambdaXiPlusFeeddownMatrix", {HistType::kTH2F, {{antilambdaPtAxis}, {antilambdaPtAxis}}});
      rFeeddownMatrices.add("hAntiLambdaAntiXiZeroFeeddownMatrix", "hAntiLambdaAntiXiZeroFeeddownMatrix", {HistType::kTH2F, {{antilambdaPtAxis}, {antilambdaPtAxis}}});
      rFeeddownMatrices.add("hAntiLambdaAntiOmegaFeeddownMatrix", "hAntiLambdaAntiOmegaPlusFeeddownMatrix", {HistType::kTH2F, {{antilambdaPtAxis}, {antilambdaPtAxis}}});
    }

    // Particle Level Corrections
    rMCCorrections.add("hK0ShSplitDenominatorPtSpectrum", "hK0ShSplitDenominatorPtSpectrum", {HistType::kTH1D, {k0ShortPtAxis}});
    rMCCorrections.add("hLambdaSplitDenominatorPtSpectrum", "hLambdaSplitDenominatorPtSpectrum", {HistType::kTH1D, {lambdaPtAxis}});
    rMCCorrections.add("hAntilambdaSplitDenominatorPtSpectrum", "hAntilambdaSplitDenominatorPtSpectrum", {HistType::kTH1F, {{antilambdaPtAxis}}});
    rMCCorrections.add("hK0ShSplitNumenatorPtSpectrum", "hK0ShSplitNumenatorPtSpectrum", {HistType::kTH1D, {k0ShortPtAxis}});
    rMCCorrections.add("hLambdaSplitNumenatorPtSpectrum", "hLambdaSplitNumenatorPtSpectrum", {HistType::kTH1D, {lambdaPtAxis}});
    rMCCorrections.add("hAntilambdaSplitNumenatorPtSpectrum", "hAntilambdaSplitNumenatorPtSpectrum", {HistType::kTH1F, {{antilambdaPtAxis}}});
    rMCCorrections.add("hK0ShBeforeEventSelectionPtSpectrum", "hK0ShBeforeEventSelectionPtSpectrum", {HistType::kTH1D, {k0ShortPtAxis}});
    rMCCorrections.add("hLambdaBeforeEventSelectionPtSpectrum", "hLambdaBeforeEventSelectionPtSpectrum", {HistType::kTH1D, {lambdaPtAxis}});
    rMCCorrections.add("hAntilambdaBeforeEventSelectionPtSpectrum", "hAntilambdaBeforeEventSelectionPtSpectrum", {HistType::kTH1F, {{antilambdaPtAxis}}});
    rMCCorrections.add("hK0ShAfterEventSelectionPtSpectrum", "hK0ShAfterEventSelectionPtSpectrum", {HistType::kTH1D, {k0ShortPtAxis}});
    rMCCorrections.add("hLambdaAfterEventSelectionPtSpectrum", "hLambdaAfterEventSelectionPtSpectrum", {HistType::kTH1D, {lambdaPtAxis}});
    rMCCorrections.add("hAntilambdaAfterEventSelectionPtSpectrum", "hAntilambdaAfterEventSelectionPtSpectrum", {HistType::kTH1F, {{antilambdaPtAxis}}});

    // Event and V0s Corrections
    rMCCorrections.add("hNEvents_Corrections", "hNEvents_Corrections", {HistType::kTH1D, {{10, 0.f, 10.f}}});
    rMCCorrections.add("hNRecEvents_MC", "hNRecEvents_MC", {HistType::kTH1D, {{1, 0.f, 1.f}}});

    // Generated Level Pt Spectrums (with rapidity cut)
    rMCCorrections.add("GenParticleRapidity", "GenParticleRapidity", {HistType::kTH1F, {{nBins, -10.0f, 10.0f}}});
    rMCCorrections.add("hK0ShGeneratedPtSpectrum", "hK0ShGeneratedPtSpectrum", {HistType::kTH1F, {k0ShortPtAxis}});
    rMCCorrections.add("hLambdaGeneratedPtSpectrum", "hLambdaGeneratedPtSpectrum", {HistType::kTH1F, {lambdaPtAxis}});
    rMCCorrections.add("hAntilambdaGeneratedPtSpectrum", "hAntilambdaGeneratedPtSpectrum", {HistType::kTH1F, {{antilambdaPtAxis}}});
    rMCCorrections.add("hXiMinusGeneratedPtSpectrum", "hXiMinusGeneratedPtSpectrum", {HistType::kTH1F, {lambdaPtAxis}});
    rMCCorrections.add("hXiZeroGeneratedPtSpectrum", "hXiZeroGeneratedPtSpectrum", {HistType::kTH1F, {lambdaPtAxis}});
    rMCCorrections.add("hOmegaGeneratedPtSpectrum", "hOmegaGeneratedPtSpectrum", {HistType::kTH1F, {lambdaPtAxis}});
    rMCCorrections.add("hXiPlusGeneratedPtSpectrum", "hXiPlusGeneratedPtSpectrum", {HistType::kTH1F, {antilambdaPtAxis}});
    rMCCorrections.add("hAntiXiZeroGeneratedPtSpectrum", "hAntiXiZeroGeneratedPtSpectrum", {HistType::kTH1F, {antilambdaPtAxis}});
    rMCCorrections.add("hAntiOmegaGeneratedPtSpectrum", "hAntiOmegaGeneratedPtSpectrum", {HistType::kTH1F, {antilambdaPtAxis}});
    rMCCorrections.add("hPhiGeneratedPtSpectrum", "hPhiGeneratedPtSpectrum", {HistType::kTH1F, {k0ShortPtAxis}});
  }

  // Event selection function
  template <typename TCollision>
  bool acceptEvent(TCollision const& collision)
  {
    rPtAnalysis.fill(HIST("hNEvents"), 0.5);
    rPtAnalysis.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(1, "All");
    if (!(collision.sel8() && dosel8)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 1.5);
    rPtAnalysis.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(2, "sel 8");
    if (!(collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && doNoTimeFrameBorder)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 2.5);
    rPtAnalysis.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(3, "NoTimeFrameBorder");
    if (!(collision.selection_bit(aod::evsel::kNoITSROFrameBorder) && doNoITSROFrameBorder)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 3.5);
    rPtAnalysis.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(4, "NoITSROFrameBorder");
    if (!(collision.selection_bit(aod::evsel::kIsTriggerTVX) && doIsTriggerTVX)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 4.5);
    rPtAnalysis.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(5, "IsTriggerTVX");
    if (!(std::abs(collision.posZ()) < cutZVertex && docutZVertex)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 5.5);
    rPtAnalysis.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(6, "cutZVertex");
    if (!(collision.selection_bit(aod::evsel::kIsVertexTOFmatched) && doIsVertexTOFmatched)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 6.5);
    rPtAnalysis.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(7, "IsVertexTOFmatched");
    if (!(collision.selection_bit(aod::evsel::kNoSameBunchPileup) && doNoSameBunchPileup)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 7.5);
    rPtAnalysis.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(8, "NoSameBunchPileup");
    if (!(collision.selection_bit(aod::evsel::kIsVertexITSTPC) && doIsVertexITSTPC)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 8.5);
    rPtAnalysis.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(9, "IsVertexITSTPC");
    if (!(collision.isInelGt0() && doisInelGt0)) {
      return false;
    }
    rPtAnalysis.fill(HIST("hNEvents"), 9.5);
    rPtAnalysis.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(10, "isInelGt0");
    return true;

    // Cut Plots
    rPtAnalysis.fill(HIST("hVertexZRec"), collision.posZ());
  }

  // V0 selection function
  template <typename TV0>
  bool acceptV0(TV0 const& v0)
  {
    const auto& posDaughterTrack = v0.template posTrack_as<DaughterTracks>(); // Positive Daughter track
    const auto& negDaughterTrack = v0.template negTrack_as<DaughterTracks>(); // Negative Daughter track

    rPtAnalysis.fill(HIST("hNV0s"), 0.5);
    rPtAnalysis.get<TH1>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(1, "All V0s");
    if (std::abs(v0.y()) > rapidityCut && doRapidityCut) { // V0 Rapidity Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNV0s"), 1.5);
    rPtAnalysis.get<TH1>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(2, "Rapidity");
    if ((std::abs(posDaughterTrack.eta()) > etadau && std::abs(negDaughterTrack.eta()) > etadau) && doDaughterPseudorapidityCut) { // Daughters Pseudorapidity Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNV0s"), 2.5);
    rPtAnalysis.get<TH1>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(3, "Dau Pseudorapidity");
    if ((posDaughterTrack.isITSAfterburner() || negDaughterTrack.isITSAfterburner()) && !doisITSAfterburner) { // ITS After Burner on daughter tracks
      return false;
    }
    rPtAnalysis.fill(HIST("hNV0s"), 3.5);
    rPtAnalysis.get<TH1>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(4, "ITS Afterburner");
    if (posDaughterTrack.itsNCls() <= itsMinHits && negDaughterTrack.itsNCls() <= itsMinHits && doitsMinHits) { // Minimum hits in the ITS
      return false;
      rPtAnalysis.fill(HIST("hNV0s"), 4.5);
      rPtAnalysis.get<TH1>(HIST("hNV0s"))->GetXaxis()->SetBinLabel(5, "ITS Min Hits");
      // Cut Plots
      rPtAnalysis.fill(HIST("V0Rapidity"), v0.y());
      rPtAnalysis.fill(HIST("hV0EtaDaughters"), v0.template posTrack_as<DaughterTracks>().eta());
      rPtAnalysis.fill(HIST("hV0EtaDaughters"), v0.template negTrack_as<DaughterTracks>().eta());
    }
    return true;
  }

  // K0sh selection function
  template <typename TV0>
  bool acceptK0sh(TV0 const& v0)
  {
    const auto& posDaughterTrack = v0.template posTrack_as<DaughterTracks>(); // Positive Daughter track
    const auto& negDaughterTrack = v0.template negTrack_as<DaughterTracks>(); // Negative Daughter track

    rPtAnalysis.fill(HIST("hNK0sh"), 0.5);
    rPtAnalysis.get<TH1>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(1, "All");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 0.5, v0.mK0Short());
    if ((std::abs(posDaughterTrack.tpcNSigmaPi()) > nSigmaTPCPion && std::abs(negDaughterTrack.tpcNSigmaPi()) > nSigmaTPCPion) && doK0shTPCPID) { // TPC PID for two pions
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 1.5);
    rPtAnalysis.get<TH1>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(2, "TPC_PID");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 1.5, v0.mK0Short());
    if (std::abs(v0.mLambda() - o2::constants::physics::MassLambda0) < compv0masscut && std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda0) < compv0masscut && doK0shcomptmasscut) { // Kzero competitive v0 mass cut (cut out Lambdas and Anti-Lambdas)
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 2.5);
    rPtAnalysis.get<TH1>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(3, "Compt_Mass");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 2.5, v0.mK0Short());
    if (v0.v0radius() > kaonshmaxct && doK0shMaxct) { // K0sh max ct
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 3.5);
    rPtAnalysis.get<TH1>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(4, "Max_ct");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 3.5, v0.mK0Short());
    if (v0.qtarm() < (k0shparamArmenterosCut * std::abs(v0.alpha())) && doK0shArmenterosCut) { // K0sh Armenteros Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 4.5);
    rPtAnalysis.get<TH1>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(5, "Armenteros");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 4.5, v0.mK0Short());
    if (v0.v0cosPA() < kaonshSettingcosPA && doK0shcosPACut) { // K0sh cosPA Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 5.5);
    rPtAnalysis.get<TH1>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(6, "cosPA");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 5.5, v0.mK0Short());
    if (v0.dcaV0daughters() > kaonshSettingdcav0dau && doK0shDCAdauCut) { // K0sh DCAdaughters Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 6.5);
    rPtAnalysis.get<TH1>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(7, "DCAdau");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 6.5, v0.mK0Short());
    if (v0.v0radius() < kaonshSettingradius && doK0shv0radiusCut) { // K0sh v0radius Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 7.5);
    rPtAnalysis.get<TH1>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(8, "v0radius");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 7.5, v0.mK0Short());
    if (std::abs(v0.dcapostopv()) < kaonshSettingdcapostopv && doK0shdcaposdautopv) { // K0sh DCAPosDaughterToPV Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 8.5);
    rPtAnalysis.get<TH1>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(9, "DCAPosDautoPV");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 8.5, v0.mK0Short());
    if (std::abs(v0.dcanegtopv()) < kaonshSettingdcanegtopv && doK0shdcanegdautopv) { // K0sh DCANegDaughterToPV Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNK0sh"), 9.5);
    rPtAnalysis.get<TH1>(HIST("hNK0sh"))->GetXaxis()->SetBinLabel(10, "DCANegDautoPV");
    rPtAnalysis.fill(HIST("hMassK0ShortvsCuts"), 9.5, v0.mK0Short());

    // Cut Plots
    rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlotK0sh"), v0.alpha(), v0.qtarm());
    rPtAnalysis.fill(HIST("hNSigmaPosPionFromK0s"), posDaughterTrack.tpcNSigmaPi(), posDaughterTrack.tpcInnerParam());
    rPtAnalysis.fill(HIST("hNSigmaNegPionFromK0s"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
    rPtAnalysis.fill(HIST("hK0shcosPA"), v0.v0cosPA());
    rPtAnalysis.fill(HIST("hK0shV0radius"), v0.v0radius());
    rPtAnalysis.fill(HIST("hK0shDCAV0Daughters"), v0.dcaV0daughters());
    rPtAnalysis.fill(HIST("hK0shDCAPosDaughter"), v0.dcapostopv());
    rPtAnalysis.fill(HIST("hK0shDCANegDaughter"), v0.dcanegtopv());
    return true;
  }

  // Lambda selection function
  template <typename TV0>
  bool acceptLambda(TV0 const& v0)
  {
    const auto& posDaughterTrack = v0.template posTrack_as<DaughterTracks>(); // Positive Daughter track
    const auto& negDaughterTrack = v0.template negTrack_as<DaughterTracks>(); // Negative Daughter track

    rPtAnalysis.fill(HIST("hNLambda"), 0.5);
    rPtAnalysis.get<TH1>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(1, "All");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 0.5, v0.mLambda());
    if (std::abs(posDaughterTrack.tpcNSigmaPr()) > nSigmaTPCProton && std::abs(negDaughterTrack.tpcNSigmaPi()) > nSigmaTPCPion && doLambdaTPCPID) { // TPC PID on daughter pion and proton for Lambda
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 1.5);
    rPtAnalysis.get<TH1>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(2, "TPC_PID");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 1.5, v0.mLambda());
    if (std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < compv0masscut && doLambdacomptmasscut) { // Lambda competitive v0 mass cut (cut out Kaons)
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 2.5);
    rPtAnalysis.get<TH1>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(3, "Compt_Mass");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 2.5, v0.mLambda());
    if (v0.v0radius() > lambdamaxct && doLambdaMaxct) { // Lambda max ct
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 3.5);
    rPtAnalysis.get<TH1>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(4, "Max_ct");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 3.5, v0.mLambda());
    if (v0.qtarm() < (lambdaparamArmenterosCut * std::abs(v0.alpha())) && doLambdaArmenterosCut) { // Lambda Armenteros Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 4.5);
    rPtAnalysis.get<TH1>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(5, "Armenteros");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 4.5, v0.mLambda());
    if (v0.v0cosPA() < lambdaSettingcosPA && doLambdacosPACut) { // Lambda cosPA Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 5.5);
    rPtAnalysis.get<TH1>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(6, "cosPA");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 5.5, v0.mLambda());
    if (v0.dcaV0daughters() > lambdaSettingdcav0dau && doLambdaDCAdauCut) { // Lambda DCAdaughters Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 6.5);
    rPtAnalysis.get<TH1>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(7, "DCAdau");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 6.5, v0.mLambda());
    if (v0.v0radius() < lambdaSettingradius && doLambdav0radiusCut) { // Lambda v0radius Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 7.5);
    rPtAnalysis.get<TH1>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(8, "v0radius");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 7.5, v0.mLambda());
    if (std::abs(v0.dcapostopv()) < lambdaSettingdcapostopv && doLambdadcaposdautopv) { // Lambda DCAPosDaughterToPV Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 8.5);
    rPtAnalysis.get<TH1>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(9, "DCAPosDautoPV");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 8.5, v0.mLambda());
    if (std::abs(v0.dcanegtopv()) < lambdaSettingdcanegtopv && doLambdadcanegdautopv) { // Lambda DCANegDaughterToPV Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNLambda"), 9.5);
    rPtAnalysis.get<TH1>(HIST("hNLambda"))->GetXaxis()->SetBinLabel(10, "DCANegDautoPV");
    rPtAnalysis.fill(HIST("hMassLambdavsCuts"), 9.5, v0.mLambda());

    // Cut Plots
    rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlotLambda"), v0.alpha(), v0.qtarm());
    rPtAnalysis.fill(HIST("hNSigmaPosProtonFromLambdas"), posDaughterTrack.tpcNSigmaPr(), posDaughterTrack.tpcInnerParam());
    rPtAnalysis.fill(HIST("hNSigmaNegPionFromLambdas"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
    rPtAnalysis.fill(HIST("hLambdacosPA"), v0.v0cosPA());
    rPtAnalysis.fill(HIST("hLambdaV0radius"), v0.v0radius());
    rPtAnalysis.fill(HIST("hLambdaDCAV0Daughters"), v0.dcaV0daughters());
    rPtAnalysis.fill(HIST("hLambdaDCAPosDaughter"), v0.dcapostopv());
    rPtAnalysis.fill(HIST("hLambdaDCANegDaughter"), v0.dcanegtopv());
    return true;
  }

  // Antilambda selection function
  template <typename TV0>
  bool acceptAntilambda(TV0 const& v0)
  {
    const auto& posDaughterTrack = v0.template posTrack_as<DaughterTracks>(); // Positive Daughter track
    const auto& negDaughterTrack = v0.template negTrack_as<DaughterTracks>(); // Negative Daughter track

    rPtAnalysis.fill(HIST("hNAntilambda"), 0.5);
    rPtAnalysis.get<TH1>(HIST("hNAntilambda"))->GetXaxis()->SetBinLabel(1, "All");
    rPtAnalysis.fill(HIST("hMassAntilambdavsCuts"), 0.5, v0.mAntiLambda());
    if (std::abs(negDaughterTrack.tpcNSigmaPr()) > nSigmaTPCProton && std::abs(posDaughterTrack.tpcNSigmaPi()) > nSigmaTPCPion) { // TPC PID on daughter pion and proton for AntiLambda
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntilambda"), 1.5);
    rPtAnalysis.get<TH1>(HIST("hNAntilambda"))->GetXaxis()->SetBinLabel(2, "TPC_PID");
    rPtAnalysis.fill(HIST("hMassAntilambdavsCuts"), 1.5, v0.mAntiLambda());
    if (std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short) < compv0masscut && doAntilambdacomptmasscut) { // Antilambda competitive v0 mass cut (cut out Kaons)
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntilambda"), 2.5);
    rPtAnalysis.get<TH1>(HIST("hNAntilambda"))->GetXaxis()->SetBinLabel(3, "Compt_Mass");
    rPtAnalysis.fill(HIST("hMassAntilambdavsCuts"), 2.5, v0.mAntiLambda());
    if (v0.v0radius() > antilambdamaxct && doAntilambdaMaxct) { // Antilambda max ct
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntilambda"), 3.5);
    rPtAnalysis.get<TH1>(HIST("hNAntilambda"))->GetXaxis()->SetBinLabel(4, "Max_ct");
    rPtAnalysis.fill(HIST("hMassAntilambdavsCuts"), 3.5, v0.mAntiLambda());
    if (v0.qtarm() < (antilambdaparamArmenterosCut * std::abs(v0.alpha())) && doAntilambdaArmenterosCut) { // Antilambda Armenteros Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntilambda"), 4.5);
    rPtAnalysis.get<TH1>(HIST("hNAntilambda"))->GetXaxis()->SetBinLabel(5, "Armenteros");
    rPtAnalysis.fill(HIST("hMassAntilambdavsCuts"), 4.5, v0.mAntiLambda());
    if (v0.v0cosPA() < antilambdaSettingcosPA && doAntilambdacosPACut) { // Antilambda cosPA Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntilambda"), 5.5);
    rPtAnalysis.get<TH1>(HIST("hNAntilambda"))->GetXaxis()->SetBinLabel(6, "cosPA");
    rPtAnalysis.fill(HIST("hMassAntilambdavsCuts"), 5.5, v0.mAntiLambda());
    if (v0.dcaV0daughters() > antilambdaSettingdcav0dau && doAntilambdaDCAdauCut) { // Antilambda DCAdaughters Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntilambda"), 6.5);
    rPtAnalysis.get<TH1>(HIST("hNAntilambda"))->GetXaxis()->SetBinLabel(7, "DCAdau");
    rPtAnalysis.fill(HIST("hMassAntilambdavsCuts"), 6.5, v0.mAntiLambda());
    if (v0.v0radius() < antilambdaSettingradius && doAntilambdav0radiusCut) { // Antilambda v0radius Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntilambda"), 7.5);
    rPtAnalysis.get<TH1>(HIST("hNAntilambda"))->GetXaxis()->SetBinLabel(8, "v0radius");
    rPtAnalysis.fill(HIST("hMassAntilambdavsCuts"), 7.5, v0.mAntiLambda());
    if (std::abs(v0.dcapostopv()) < antilambdaSettingdcapostopv && doAntilambdadcaposdautopv) { // Antilambda DCAPosDaughterToPV Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntilambda"), 8.5);
    rPtAnalysis.get<TH1>(HIST("hNAntilambda"))->GetXaxis()->SetBinLabel(9, "DCAPosDautoPV");
    rPtAnalysis.fill(HIST("hMassAntilambdavsCuts"), 8.5, v0.mAntiLambda());
    if (std::abs(v0.dcanegtopv()) < antilambdaSettingdcanegtopv && doAntilambdadcanegdautopv) { // Antilambda DCANegDaughterToPV Topological Cut
      return false;
    }
    rPtAnalysis.fill(HIST("hNAntilambda"), 9.5);
    rPtAnalysis.get<TH1>(HIST("hNAntilambda"))->GetXaxis()->SetBinLabel(10, "DCANegDautoPV");
    rPtAnalysis.fill(HIST("hMassAntilambdavsCuts"), 9.5, v0.mAntiLambda());

    // Cut plots
    rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlotAntilambda"), v0.alpha(), v0.qtarm());
    rPtAnalysis.fill(HIST("hNSigmaPosPionFromAntilambdas"), posDaughterTrack.tpcNSigmaPr(), posDaughterTrack.tpcInnerParam());
    rPtAnalysis.fill(HIST("hNSigmaNegProtonFromAntilambdas"), negDaughterTrack.tpcNSigmaPi(), negDaughterTrack.tpcInnerParam());
    rPtAnalysis.fill(HIST("hAntilambdacosPA"), v0.v0cosPA());
    rPtAnalysis.fill(HIST("hAntilambdaV0radius"), v0.v0radius());
    rPtAnalysis.fill(HIST("hAntilambdaDCAV0Daughters"), v0.dcaV0daughters());
    rPtAnalysis.fill(HIST("hAntilambdaDCAPosDaughter"), v0.dcapostopv());
    rPtAnalysis.fill(HIST("hAntilambdaDCANegDaughter"), v0.dcanegtopv());
    return true;
  }

  // Defining the type of the daughter tracks
  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi, aod::pidTPCPr>;

  void genMCProcess(
    aod::McCollisions::iterator const& /*mcCollisions*/,
    soa::SmallGroups<soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, /*aod::McCentFT0Ms,*/ aod::PVMults>> const& collisions,
    aod::McParticles const& mcParticles)
  {
    // Event Efficiency, Event Split and V0 Signal Loss Corrections
    rMCCorrections.fill(HIST("hNEvents_Corrections"), 0.5); // Event Efficiency Denominator

    // Particles (of interest) Generated Pt Spectrum and Signal Loss Denominator Loop
    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) < rapidityCut) {
        if (mcParticle.isPhysicalPrimary()) {
          rMCCorrections.fill(HIST("GenParticleRapidity"), mcParticle.y());
          if (mcParticle.pdgCode() == kK0Short) // kzero matched
          {
            rMCCorrections.fill(HIST("hK0ShGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == kLambda0) // lambda matched
          {
            rMCCorrections.fill(HIST("hLambdaGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == kLambda0Bar) // antilambda matched
          {
            rMCCorrections.fill(HIST("hAntilambdaGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == kXiMinus) // Xi Minus matched
          {
            rMCCorrections.fill(HIST("hXiMinusGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == kXi0) // Xi Zero matched
          {
            rMCCorrections.fill(HIST("hXiZeroGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == kOmegaMinus) // Omega matched
          {
            rMCCorrections.fill(HIST("hOmegaGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == kXiPlusBar) // Xi Plus matched
          {
            rMCCorrections.fill(HIST("hXiPlusGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == -kXi0) // Anti-Xi Zero matched
          {
            rMCCorrections.fill(HIST("hAntiXiZeroGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == kOmegaPlusBar) // Anti-Omega matched
          {
            rMCCorrections.fill(HIST("hAntiOmegaGeneratedPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == kPhi) // Phi
          {
            rMCCorrections.fill(HIST("hPhiGeneratedPtSpectrum"), mcParticle.pt());
          }
        }
      }
    }
    // Signal Loss Numenator Loop
    for (const auto& collision : collisions) {
      rMCCorrections.fill(HIST("hNEvents_Corrections"), 1.5); // Number of Events Reconsctructed
      if (!acceptEvent(collision)) {                          // Event Selection
        return;
      }
      rMCCorrections.fill(HIST("hNEvents_Corrections"), 2.5); // Event Split Denomimator and Event Efficiency Numenator
      for (const auto& mcParticle : mcParticles) {
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }
        if (std::abs(mcParticle.y()) > rapidityCut) {
          continue;
        }
        if (mcParticle.pdgCode() == kK0Short) // kzero matched
        {
          rMCCorrections.fill(HIST("hK0ShAfterEventSelectionPtSpectrum"), mcParticle.pt());
        }
        if (mcParticle.pdgCode() == kLambda0) // lambda matched
        {
          rMCCorrections.fill(HIST("hLambdaAfterEventSelectionPtSpectrum"), mcParticle.pt());
        }
        if (mcParticle.pdgCode() == kLambda0Bar) // antilambda matched
        {
          rMCCorrections.fill(HIST("hAntilambdaAfterEventSelectionPtSpectrum"), mcParticle.pt());
        }
      }
    }
    // End of Signal Loss Numenator Loop
  }
  // This is the Process for the MC reconstructed Data
  void recMCProcess(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::PVMults>::iterator const& collision,
                    soa::Join<aod::V0Datas, aod::McV0Labels> const& V0s,
                    DaughterTracks const&, // no need to define a variable for tracks, if we don't access them directly
                    aod::McParticles const& mcParticles)
  {
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
    if (!acceptEvent(collision)) { // Event Selection
      return;
    }
    rMCCorrections.fill(HIST("hNRecEvents_MC"), 0.5); // Event Split Numenator

    // v0 Signal Splitting Numenator Start
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary()) {
        if (std::abs(mcParticle.y()) < rapidityCut) {
          if (mcParticle.pdgCode() == kK0Short) { // kzero matched
            rMCCorrections.fill(HIST("hK0ShSplitNumenatorPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == kLambda0) { // lambda matched
            rMCCorrections.fill(HIST("hLambdaSplitNumenatorPtSpectrum"), mcParticle.pt());
          }
          if (mcParticle.pdgCode() == kLambda0Bar) { // antilambda matched
            rMCCorrections.fill(HIST("hAntilambdaSplitNumenatorPtSpectrum"), mcParticle.pt());
          }
        }
      }
    }
    // V0 Signal Splitting Numenator End

    for (const auto& v0 : V0s) {
      // Checking that the V0 is a true K0s/Lambdas/Antilambdas and then filling the parameter histograms and the invariant mass plots for different cuts (which are taken from namespace)
      if (v0.has_mcParticle()) {
        auto v0mcParticle = v0.mcParticle();

        // signal splitting demoninator
        if (v0mcParticle.isPhysicalPrimary()) {
          if (v0mcParticle.pdgCode() == kK0Short) { // kzero matched
            rMCCorrections.fill(HIST("hK0ShSplitDenominatorPtSpectrum"), v0mcParticle.pt());
          }
          if (v0mcParticle.pdgCode() == kLambda0) { // lambda matched
            rMCCorrections.fill(HIST("hLambdaSplitDenominatorPtSpectrum"), v0mcParticle.pt());
          }
          if (v0mcParticle.pdgCode() == kLambda0Bar) { // antilambda matched
            rMCCorrections.fill(HIST("hAntilambdaSplitDenominatorPtSpectrum"), v0mcParticle.pt());
          }
        }
        // signal splitting demoninator end

        if (!acceptV0(v0)) { // V0 Selections
          continue;
        }
        rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlot"), v0.alpha(), v0.qtarm());
        // kzero analysis
        if (kzeroAnalysis == true) {
          if (v0mcParticle.pdgCode() == kK0Short) { // kzero matched
            if (!acceptK0sh(v0)) {                  // K0sh Selection
              continue;
            }

            if (v0mcParticle.isPhysicalPrimary()) {
              for (int i = 0; i < nmaxHistograms; i++) {
                if (kaonptedgevalues[i] <= v0.pt() && v0.pt() < kaonptedgevalues[i + 1]) { // finding v0s with pt within the range of our bin edges
                  pthistos::kaonPt[i]->Fill(v0.mK0Short());                                // filling the k0s namespace histograms
                }
              }
            }
            if (!v0mcParticle.isPhysicalPrimary()) {
              auto v0mothers = v0mcParticle.mothers_as<aod::McParticles>(); // Get mothers
              if (!v0mothers.empty()) {
                auto& v0mcParticleMother = v0mothers.front(); // First mother
                rFeeddownMatrices.fill(HIST("hK0shFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                if (v0mcParticleMother.pdgCode() == kPhi) { // Phi Mother Matched
                  rFeeddownMatrices.fill(HIST("hK0shPhiFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                }
              }
            }
          }
        }
        // lambda analysis
        if (lambdaAnalysis == true) {
          if (v0mcParticle.pdgCode() == kLambda0) { // lambda matched

            if (!acceptLambda(v0)) { // Lambda Selections
              continue;
            }

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
                rFeeddownMatrices.fill(HIST("hLambdaFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                if (v0mcParticleMother.pdgCode() == kXiMinus) { // Xi Minus Mother Matched
                  rFeeddownMatrices.fill(HIST("hLambdaXiMinusFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                }
                if (v0mcParticleMother.pdgCode() == kXi0) { // Xi Zero Mother Matched
                  rFeeddownMatrices.fill(HIST("hLambdaXiZeroFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                }
                if (v0mcParticleMother.pdgCode() == kOmegaMinus) { // Omega Mother Matched
                  rFeeddownMatrices.fill(HIST("hLambdaOmegaFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                }
              }
            }
          }
        }
        // antilambda analysis
        if (antiLambdaAnalysis == true) {
          if (v0mcParticle.pdgCode() == kLambda0Bar) { // antilambda matched

            if (!acceptAntilambda(v0)) { // Antilambda Selections
              continue;
            }

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
                rFeeddownMatrices.fill(HIST("hAntiLambdaFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                if (v0mcParticleMother.pdgCode() == kXiPlusBar) { // Xi Plus Mother Matched
                  rFeeddownMatrices.fill(HIST("hAntiLambdaXiPlusFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                }
                if (v0mcParticleMother.pdgCode() == -kXi0) { // Anti-Xi Zero Mother Matched
                  rFeeddownMatrices.fill(HIST("hAntiLambdaAntiXiZeroFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                }
                if (v0mcParticleMother.pdgCode() == kOmegaPlusBar) { // Anti-Omega (minus) Mother Matched
                  rFeeddownMatrices.fill(HIST("hAntiLambdaAntiOmegaFeeddownMatrix"), v0mcParticle.pt(), v0mcParticleMother.pt());
                }
              }
            }
          }
        }
      }
    }
  }
  // This is the process for Real Data
  void dataProcess(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults>::iterator const& collision,
                   aod::V0Datas const& V0s,
                   DaughterTracks const&)
  {
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
    if (!acceptEvent(collision)) { // Event Selection
      return;
    }
    rPtAnalysis.fill(HIST("hNRecEvents_Data"), 1.0); // Number of Reconstructed Events

    for (const auto& v0 : V0s) {
      // Checking that the V0 is a true K0s/Lambdas/Antilambdas and then filling the parameter histograms and the invariant mass plots for different cuts (which are taken from namespace)
      if (!acceptV0(v0)) { // V0 Selection
        continue;
      }
      rPtAnalysis.fill(HIST("hArmenterosPodolanskiPlot"), v0.alpha(), v0.qtarm());
      // kzero analysis
      if (kzeroAnalysis == true) {
        if (!acceptK0sh(v0)) { // K0sh Selection
          continue;
        }
        for (int i = 0; i < nmaxHistograms; i++) {
          if (kaonptedgevalues[i] <= v0.pt() && v0.pt() < kaonptedgevalues[i + 1]) { // finding v0s with pt within the range of our bin edges
            pthistos::kaonPt[i]->Fill(v0.mK0Short());                                // filling the k0s namespace histograms
          }
        }
      }
      // lambda analysis
      if (lambdaAnalysis == true) {
        if (!acceptLambda(v0)) { // Lambda Selection
          continue;
        }
        for (int i = 0; i < nmaxHistograms; i++) {
          if (lambdaptedgevalues[i] <= v0.pt() && v0.pt() < lambdaptedgevalues[i + 1]) {
            pthistos::lambdaPt[i]->Fill(v0.mLambda());
          }
        }
      }
      // anti-lambda analysis
      if (antiLambdaAnalysis == true) {
        if (!acceptAntilambda(v0)) { // Antilambda Selection
          continue;
        }
        for (int i = 0; i < nmaxHistograms; i++) {
          if (lambdaptedgevalues[i] <= v0.pt() && v0.pt() < lambdaptedgevalues[i + 1]) {
            pthistos::antilambdaPt[i]->Fill(v0.mAntiLambda());
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
