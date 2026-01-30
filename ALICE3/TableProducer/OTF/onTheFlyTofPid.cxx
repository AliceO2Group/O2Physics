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
/// \file onTheFlyTofPid.cxx
///
/// \brief This task goes straight from a combination of track table and mcParticles
///        and a custom TOF configuration to a table of TOF NSigmas for the particles
///        being analysed. It currently contemplates 5 particle types:
///        electrons, pions, kaons, protons and muons
///
///        More particles could be added but would have to be added to the LUT
///        being used in the onTheFly tracker task.
///
/// \author David Dobrigkeit Chinellato, UNICAMP
/// \author Nicola Nicassio, University and INFN Bari
/// \since  May 22, 2024
///

#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "ALICE3/Core/FastTracker.h"
#include "ALICE3/Core/TrackUtilities.h"
#include "ALICE3/DataModel/OTFCollision.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/GeomConstants.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <CommonUtils/NameConf.h>
#include <DataFormatsCalibration/MeanVertexObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/HelixHelper.h>

#include <TEfficiency.h>
#include <THashList.h>
#include <TPDGCode.h>
#include <TRandom3.h>

#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

std::array<std::shared_ptr<TH2>, 9> h2dInnerTimeResTrack;
std::array<std::shared_ptr<TH2>, 9> h2dInnerTimeResTotal;
std::array<std::shared_ptr<TH2>, 9> h2dOuterTimeResTrack;
std::array<std::shared_ptr<TH2>, 9> h2dOuterTimeResTotal;
std::array<std::array<std::shared_ptr<TH2>, 9>, 9> h2dInnerNsigmaTrue;
std::array<std::array<std::shared_ptr<TH2>, 9>, 9> h2dOuterNsigmaTrue;
std::array<std::array<std::shared_ptr<TH2>, 9>, 9> h2dInnerDeltaTrue;
std::array<std::array<std::shared_ptr<TH2>, 9>, 9> h2dOuterDeltaTrue;

struct OnTheFlyTofPid {
  Produces<aod::UpgradeTofMC> upgradeTofMC;
  Produces<aod::UpgradeTof> upgradeTof;
  Produces<aod::UpgradeTofExpectedTime> upgradeTofExpectedTime;

  // necessary for particle charges
  Service<o2::framework::O2DatabasePDG> pdg;
  // Necessary for LUTs
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // these are the settings governing the TOF layers to be used
  // note that there are two layers foreseen for now: inner and outer TOF
  // more could be added (especially a disk TOF at a certain z?)
  // in the evolution of this effort
  struct : ConfigurableGroup {
    Configurable<bool> considerEventTime{"considerEventTime", true, "flag to consider event time"};
    Configurable<float> innerTOFRadius{"innerTOFRadius", 20, "barrel inner TOF radius (cm)"};
    Configurable<float> outerTOFRadius{"outerTOFRadius", 80, "barrel outer TOF radius (cm)"};
    Configurable<float> innerTOFLength{"innerTOFLength", 124, "barrel inner TOF length (cm)"};
    Configurable<float> outerTOFLength{"outerTOFLength", 250, "barrel outer TOF length (cm)"};
    Configurable<std::array<float, 2>> innerTOFPixelDimension{"innerTOFPixelDimension", {0.1, 0.1}, "barrel inner TOF pixel dimension in Z and RPhi (cm)"};
    Configurable<float> innerTOFFractionOfInactiveArea{"innerTOFFractionOfInactiveArea", 0.f, "barrel inner TOF fraction of inactive area"};
    Configurable<std::array<float, 2>> outerTOFPixelDimension{"outerTOFPixelDimension", {0.1, 0.1}, "barrel outer TOF pixel dimension in Z and RPhi (cm)"};
    Configurable<float> outerTOFFractionOfInactiveArea{"outerTOFFractionOfInactiveArea", 0.f, "barrel outer TOF fraction of inactive area"};
    Configurable<float> innerTOFTimeReso{"innerTOFTimeReso", 20, "barrel inner TOF time error (ps)"};
    Configurable<float> outerTOFTimeReso{"outerTOFTimeReso", 20, "barrel outer TOF time error (ps)"};
    Configurable<int> nStepsLIntegrator{"nStepsLIntegrator", 200, "number of steps in length integrator"};
    Configurable<float> multiplicityEtaRange{"multiplicityEtaRange", 0.800000012, "eta range to compute the multiplicity"};
    Configurable<bool> flagIncludeTrackTimeRes{"flagIncludeTrackTimeRes", true, "flag to include or exclude track time resolution"};
    Configurable<bool> flagTOFLoadDelphesLUTs{"flagTOFLoadDelphesLUTs", false, "flag to load Delphes LUTs for tracking correction (use recoTrack parameters if false)"};
  } simConfig;

  struct : ConfigurableGroup {
    Configurable<bool> doQAplots{"doQAplots", true, "do basic velocity plot qa"};
    Configurable<bool> doSeparationVsPt{"doSeparationVsPt", true, "Produce plots vs pt or p"};
    Configurable<int> nBinsBeta{"nBinsBeta", 2200, "number of bins in beta"};
    Configurable<int> nBinsP{"nBinsP", 80, "number of bins in momentum/pT depending on doSeparationVsPt"};
    Configurable<int> nBinsTrackLengthInner{"nBinsTrackLengthInner", 300, "number of bins in track length"};
    Configurable<int> nBinsTrackLengthOuter{"nBinsTrackLengthOuter", 300, "number of bins in track length"};
    Configurable<int> nBinsTrackDeltaLength{"nBinsTrackDeltaLength", 100, "number of bins in delta track length"};
    Configurable<int> nBinsNsigmaCorrectSpecies{"nBinsNsigmaCorrectSpecies", 200, "number of bins in Nsigma plot (correct speies)"};
    Configurable<int> nBinsNsigmaWrongSpecies{"nBinsNsigmaWrongSpecies", 200, "number of bins in Nsigma plot (wrong species)"};
    Configurable<float> minNsigmaRange{"minNsigmaRange", -10, "lower limit for the nsigma axis"};
    Configurable<float> maxNsigmaRange{"maxNsigmaRange", +10, "upper limit for the nsigma axis"};
    Configurable<int> nBinsDeltaCorrectSpecies{"nBinsDeltaCorrectSpecies", 200, "number of bins in Delta plot (correct speies)"};
    Configurable<int> nBinsDeltaWrongSpecies{"nBinsDeltaWrongSpecies", 200, "number of bins in Delta plot (wrong species)"};
    Configurable<float> minDeltaRange{"minDeltaRange", -100, "lower limit for the nsigma axis"};
    Configurable<float> maxDeltaRange{"maxDeltaRange", +100, "upper limit for the nsigma axis"};
    Configurable<int> nBinsTimeRes{"nBinsTimeRes", 400, "number of bins plots time resolution"};
    Configurable<int> nBinsRelativeEtaPt{"nBinsRelativeEtaPt", 400, "number of bins plots pt and eta relative errors"};
    Configurable<int> nBinsEta{"nBinsEta", 400, "number of bins plot relative eta error"};
    Configurable<int> nBinsMult{"nBinsMult", 200, "number of bins in multiplicity"};
    Configurable<float> maxMultRange{"maxMultRange", 1000.f, "upper limit in multiplicity plots"};
  } plotsConfig;
  Configurable<bool> cleanLutWhenLoaded{"cleanLutWhenLoaded", true, "clean LUTs after being loaded to save disk space"};

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Track smearer array, one per geometry
  std::vector<std::unique_ptr<o2::delphes::DelphesO2TrackSmearer>> mSmearer;

  // needed: random number generator for smearing
  TRandom3 pRandomNumberGenerator;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<THashList> listEfficiency{"efficiency"};
  static constexpr int kParticles = 9;

  // Configuration defined at init time
  o2::fastsim::GeometryContainer mGeoContainer;
  float mMagneticField = 0.0f;
  void init(o2::framework::InitContext& initContext)
  {
    mGeoContainer.init(initContext);

    const int nGeometries = mGeoContainer.getNumberOfConfigurations();
    mMagneticField = mGeoContainer.getFloatValue(0, "global", "magneticfield");

    pRandomNumberGenerator.SetSeed(0); // fully randomize

    for (int icfg = 0; icfg < nGeometries; ++icfg) {
      const std::string histPath = "Configuration_" + std::to_string(icfg) + "/";
      mSmearer.emplace_back(std::make_unique<o2::delphes::DelphesO2TrackSmearer>());
      mSmearer[icfg]->setCleanupDownloadedFile(cleanLutWhenLoaded.value);
      mSmearer[icfg]->setCcdbManager(ccdb.operator->());
      mSmearer[icfg]->setDownloadPath("TOFPID");
      std::map<std::string, std::string> globalConfiguration = mGeoContainer.getConfiguration(icfg, "global");
      for (const auto& entry : globalConfiguration) {
        int pdg = 0;
        if (entry.first.find("lut") != 0) {
          continue;
        }
        if (entry.first.find("lutEl") != std::string::npos) {
          pdg = kElectron;
        } else if (entry.first.find("lutMu") != std::string::npos) {
          pdg = kMuonMinus;
        } else if (entry.first.find("lutPi") != std::string::npos) {
          pdg = kPiPlus;
        } else if (entry.first.find("lutKa") != std::string::npos) {
          pdg = kKPlus;
        } else if (entry.first.find("lutPr") != std::string::npos) {
          pdg = kProton;
        } else if (entry.first.find("lutDe") != std::string::npos) {
          pdg = o2::constants::physics::kDeuteron;
        } else if (entry.first.find("lutTr") != std::string::npos) {
          pdg = o2::constants::physics::kTriton;
        } else if (entry.first.find("lutHe3") != std::string::npos) {
          pdg = o2::constants::physics::kHelium3;
        } else if (entry.first.find("lutAl") != std::string::npos) {
          pdg = o2::constants::physics::kAlpha;
        }

        std::string filename = entry.second;
        if (pdg == 0) {
          LOG(fatal) << "Unknown LUT entry " << entry.first << " for global configuration";
        }
        LOG(info) << "Loading LUT for pdg " << pdg << " for config " << icfg << " from provided file '" << filename << "'";
        if (filename.empty()) {
          LOG(warning) << "No LUT file passed for pdg " << pdg << ", skipping.";
        }
        // strip from leading/trailing spaces
        filename.erase(0, filename.find_first_not_of(" "));
        filename.erase(filename.find_last_not_of(" ") + 1);
        if (filename.empty()) {
          LOG(warning) << "No LUT file passed for pdg " << pdg << ", skipping.";
        }
        bool success = mSmearer[icfg]->loadTable(pdg, filename.c_str());
        if (!success) {
          LOG(fatal) << "Having issue with loading the LUT " << pdg << " " << filename;
        }
      }
    }

    if (plotsConfig.doQAplots) {
      const AxisSpec axisdNdeta{plotsConfig.nBinsMult, 0.0f, plotsConfig.maxMultRange, Form("dN/d#eta in |#eta| < %f", simConfig.multiplicityEtaRange.value)};

      histos.add("h1dNdeta", "h2dNdeta", kTH1F, {axisdNdeta});
      histos.add("h2dEventTime", "h2dEventTime", kTH2F, {{200, -1000, 1000, "computed"}, {200, -1000, 1000, "generated"}});
      histos.add("h1dEventTimegen", "h1dEventTimegen", kTH1F, {{200, -1000, 1000, "generated"}});
      histos.add("h1dEventTimerec", "h1dEventTimerec", kTH1F, {{200, -1000, 1000, "computed"}});
      histos.add("h1dEventTimedelta", "h1dEventTimedelta", kTH1F, {{200, -1000, 1000, "generated - computed"}});
      histos.add("h2dEventTimeres", "h2dEventTimeres", kTH2F, {axisdNdeta, {300, 0, 300, "resolution"}});
      listEfficiency.setObject(new THashList);
      listEfficiency->Add(new TEfficiency("effEventTime", "effEventTime", plotsConfig.nBinsMult, 0.0f, plotsConfig.maxMultRange));

      const AxisSpec axisMomentum{static_cast<int>(plotsConfig.nBinsP), 0.0f, +10.0f, "#it{p} (GeV/#it{c})"};
      const AxisSpec axisRigidity{static_cast<int>(plotsConfig.nBinsP), 0.0f, +10.0f, "#it{p} / |#it{z}| (GeV/#it{c})"};
      const AxisSpec axisMomentumSmall{static_cast<int>(plotsConfig.nBinsP), 0.0f, +1.0f, "#it{p} (GeV/#it{c})"};
      const AxisSpec axisVelocity{static_cast<int>(plotsConfig.nBinsBeta), 0.0f, +1.1f, "Measured #beta"};
      const AxisSpec axisTrackLengthInner{static_cast<int>(plotsConfig.nBinsTrackLengthInner), 0.0f, 60.0f, "Track length (cm)"};
      const AxisSpec axisTrackLengthOuter{static_cast<int>(plotsConfig.nBinsTrackLengthOuter), 0.0f, 300.0f, "Track length (cm)"};
      const AxisSpec axisTrackDeltaLength{static_cast<int>(plotsConfig.nBinsTrackDeltaLength), 0.0f, 30.0f, "Delta Track length (cm)"};
      histos.add("iTOF/h2dVelocityVsMomentumInner", "h2dVelocityVsMomentumInner", kTH2F, {axisMomentum, axisVelocity});
      histos.add("iTOF/h2dVelocityVsRigidityInner", "h2dVelocityVsRigidityInner", kTH2F, {axisRigidity, axisVelocity});
      histos.add("iTOF/h2dTrackLengthInnerVsPt", "h2dTrackLengthInnerVsPt", kTH2F, {axisMomentumSmall, axisTrackLengthInner});
      histos.add("iTOF/h2dTrackLengthInnerRecoVsPt", "h2dTrackLengthInnerRecoVsPt", kTH2F, {axisMomentumSmall, axisTrackLengthInner});
      histos.add("iTOF/h2dDeltaTrackLengthInnerVsPt", "h2dDeltaTrackLengthInnerVsPt", kTH2F, {axisMomentumSmall, axisTrackDeltaLength});
      histos.add("iTOF/h2HitMap", "h2HitMap", kTH2F, {{1000, -simConfig.innerTOFLength / 2, simConfig.innerTOFLength / 2}, {1000, 0, simConfig.innerTOFRadius * 2 * M_PI}});

      histos.add("oTOF/h2dVelocityVsMomentumOuter", "h2dVelocityVsMomentumOuter", kTH2F, {axisMomentum, axisVelocity});
      histos.add("oTOF/h2dVelocityVsRigidityOuter", "h2dVelocityVsRigidityOuter", kTH2F, {axisRigidity, axisVelocity});
      histos.add("oTOF/h2dTrackLengthOuterVsPt", "h2dTrackLengthOuterVsPt", kTH2F, {axisMomentumSmall, axisTrackLengthOuter});
      histos.add("oTOF/h2dTrackLengthOuterRecoVsPt", "h2dTrackLengthOuterRecoVsPt", kTH2F, {axisMomentumSmall, axisTrackLengthOuter});
      histos.add("oTOF/h2dDeltaTrackLengthOuterVsPt", "h2dDeltaTrackLengthOuterVsPt", kTH2F, {axisMomentumSmall, axisTrackDeltaLength});
      histos.add("oTOF/h2HitMap", "h2HitMap", kTH2F, {{1000, -simConfig.outerTOFLength / 2, simConfig.outerTOFLength / 2}, {1000, 0, simConfig.outerTOFRadius * 2 * M_PI}});

      const AxisSpec axisPt{static_cast<int>(plotsConfig.nBinsP), 0.0f, +4.0f, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec axisEta{static_cast<int>(plotsConfig.nBinsEta), -2.0f, +2.0f, "#eta"};
      const AxisSpec axisRelativePt{static_cast<int>(plotsConfig.nBinsRelativeEtaPt), 0.0f, +10.0f, "#it{#sigma_{p_{T}}} / #it{p_{T}} (%)"};
      const AxisSpec axisRelativeEta{static_cast<int>(plotsConfig.nBinsRelativeEtaPt), 0.0f, +10.0f, "#it{#sigma_{#eta}} / #it{#eta} (%)"};
      histos.add("h2dRelativePtResolution", "h2dRelativePtResolution", kTH2F, {axisPt, axisRelativePt});
      histos.add("h2dRelativeEtaResolution", "h2dRelativeEtaResolution", kTH2F, {axisEta, axisRelativeEta});

      std::string particleNames[kParticles] = {"#it{e}", "#it{#mu}", "#it{#pi}", "#it{K}", "#it{p}", "#it{d}", "#it{t}", "^{3}He", "#it{#alpha}"};
      std::string particleNames2[kParticles] = {"Elec", "Muon", "Pion", "Kaon", "Prot", "Deut", "Trit", "He3", "Al"};
      for (int iTrue = 0; iTrue < kParticles; iTrue++) {
        auto addHistogram = [&](const std::string& name, const AxisSpec& axis) {
          return histos.add<TH2>(name, "", kTH2F, {axisMomentum, axis});
        };

        const AxisSpec axisTrackTimeRes{plotsConfig.nBinsTimeRes, 0.0f, +200.0f, "Track time resolution - " + particleNames[iTrue] + " (ps)"};
        h2dInnerTimeResTrack[iTrue] = addHistogram("iTOF/res/h2dInnerTimeResTrack" + particleNames2[iTrue] + "VsP", axisTrackTimeRes);
        h2dOuterTimeResTrack[iTrue] = addHistogram("oTOF/res/h2dOuterTimeResTrack" + particleNames2[iTrue] + "VsP", axisTrackTimeRes);
        const AxisSpec axisTotalTimeRes{plotsConfig.nBinsTimeRes, 0.0f, +200.0f, "Total time resolution - " + particleNames[iTrue] + " (ps)"};
        h2dInnerTimeResTotal[iTrue] = addHistogram("iTOF/res/h2dInnerTimeResTotal" + particleNames2[iTrue] + "VsP", axisTotalTimeRes);
        h2dOuterTimeResTotal[iTrue] = addHistogram("oTOF/res/h2dOuterTimeResTotal" + particleNames2[iTrue] + "VsP", axisTotalTimeRes);
        for (int iHyp = 0; iHyp < kParticles; iHyp++) {
          std::string nameTitleInner = "h2dInnerNsigmaTrue" + particleNames2[iTrue] + "Vs" + particleNames2[iHyp] + "Hypothesis";
          std::string nameTitleOuter = "h2dOuterNsigmaTrue" + particleNames2[iTrue] + "Vs" + particleNames2[iHyp] + "Hypothesis";
          std::string nameTitleInnerDelta = "h2dInnerDeltaTrue" + particleNames2[iTrue] + "Vs" + particleNames2[iHyp] + "Hypothesis";
          std::string nameTitleOuterDelta = "h2dOuterDeltaTrue" + particleNames2[iTrue] + "Vs" + particleNames2[iHyp] + "Hypothesis";
          const AxisSpec axisX{plotsConfig.doSeparationVsPt.value ? axisPt : axisMomentum};
          const AxisSpec axisNsigmaCorrect{plotsConfig.nBinsNsigmaCorrectSpecies, plotsConfig.minNsigmaRange, plotsConfig.maxNsigmaRange, "N#sigma - True " + particleNames[iTrue] + " vs " + particleNames[iHyp] + " hypothesis"};
          const AxisSpec axisDeltaCorrect{plotsConfig.nBinsDeltaCorrectSpecies, plotsConfig.minDeltaRange, plotsConfig.maxDeltaRange, "#Delta - True " + particleNames[iTrue] + " vs " + particleNames[iHyp] + " hypothesis"};
          const AxisSpec axisNsigmaWrong{plotsConfig.nBinsNsigmaWrongSpecies, plotsConfig.minNsigmaRange, plotsConfig.maxNsigmaRange, "N#sigma -  True " + particleNames[iTrue] + " vs " + particleNames[iHyp] + " hypothesis"};
          const AxisSpec axisDeltaWrong{plotsConfig.nBinsDeltaWrongSpecies, plotsConfig.minDeltaRange, plotsConfig.maxDeltaRange, "#Delta - True " + particleNames[iTrue] + " vs " + particleNames[iHyp] + " hypothesis"};
          const AxisSpec axisNSigma{iTrue == iHyp ? axisNsigmaCorrect : axisNsigmaWrong};
          const AxisSpec axisDelta{iTrue == iHyp ? axisDeltaCorrect : axisDeltaWrong};
          h2dInnerNsigmaTrue[iTrue][iHyp] = histos.add<TH2>("iTOF/nsigma/h2dInnerNsigmaTrue" + particleNames2[iTrue] + "Vs" + particleNames2[iHyp] + "Hypothesis", "", kTH2F, {axisX, axisNSigma});
          h2dOuterNsigmaTrue[iTrue][iHyp] = histos.add<TH2>("oTOF/nsigma/h2dOuterNsigmaTrue" + particleNames2[iTrue] + "Vs" + particleNames2[iHyp] + "Hypothesis", "", kTH2F, {axisX, axisNSigma});
          h2dInnerDeltaTrue[iTrue][iHyp] = histos.add<TH2>("iTOF/delta/h2dInnerDeltaTrue" + particleNames2[iTrue] + "Vs" + particleNames2[iHyp] + "Hypothesis", "", kTH2F, {axisX, axisDelta});
          h2dOuterDeltaTrue[iTrue][iHyp] = histos.add<TH2>("oTOF/delta/h2dOuterDeltaTrue" + particleNames2[iTrue] + "Vs" + particleNames2[iHyp] + "Hypothesis", "", kTH2F, {axisX, axisDelta});
        }
      }
    }
  }

  struct TOFLayerEfficiency {
   private:
    const float layerRadius;
    const float layerLength;
    const float pixelDimensionZ;
    const float pixelDimensionRPhi;
    const float fractionInactive;
    const float magField;

    TAxis* axisZ = nullptr;
    TAxis* axisRPhi = nullptr;
    TAxis* axisInPixelZ = nullptr;
    TAxis* axisInPixelRPhi = nullptr;

    TH2F* hHitMapInPixel = nullptr;
    TH2F* hHitMapInPixelBefore = nullptr;
    TH2F* hHitMap = nullptr;

   public:
    ~TOFLayerEfficiency()
    {
      hHitMap->SaveAs(Form("/tmp/%s.png", hHitMap->GetName()));
      hHitMapInPixel->SaveAs(Form("/tmp/%s.png", hHitMapInPixel->GetName()));
      hHitMapInPixelBefore->SaveAs(Form("/tmp/%s.png", hHitMapInPixelBefore->GetName()));

      delete axisZ;
      delete axisRPhi;
      delete axisInPixelZ;
      delete axisInPixelRPhi;
      delete hHitMap;
      delete hHitMapInPixel;
      delete hHitMapInPixelBefore;
    }

    TOFLayerEfficiency(float r, float l, std::array<float, 2> pDimensions, float fIA, float m)
      : layerRadius(r), layerLength(l), pixelDimensionZ(pDimensions[0]), pixelDimensionRPhi(pDimensions[1]), fractionInactive(fIA), magField(m)
    {
      // Assuming square pixels for simplicity
      const float circumference = 2.0f * M_PI * layerRadius;
      axisZ = new TAxis(static_cast<int>(layerLength / pixelDimensionZ), -layerLength / 2, layerLength);
      axisRPhi = new TAxis(static_cast<int>(circumference / pixelDimensionRPhi), 0.f, circumference);

      const float inactiveBorderRPhi = pixelDimensionRPhi * std::sqrt(fractionInactive) / 2;
      const float inactiveBorderZ = pixelDimensionZ * std::sqrt(fractionInactive) / 2;
      const double arrayRPhi[4] = {-pixelDimensionRPhi / 2, -pixelDimensionRPhi / 2 + inactiveBorderRPhi, pixelDimensionRPhi / 2 - inactiveBorderRPhi, pixelDimensionRPhi / 2};
      for (int i = 0; i < 4; i++) {
        LOG(info) << "arrayRPhi[" << i << "] = " << arrayRPhi[i];
      }
      axisInPixelRPhi = new TAxis(3, arrayRPhi);
      const double arrayZ[4] = {-pixelDimensionZ / 2, -pixelDimensionZ / 2 + inactiveBorderZ, pixelDimensionZ / 2 - inactiveBorderZ, pixelDimensionZ / 2};
      for (int i = 0; i < 4; i++) {
        LOG(info) << "arrayZ[" << i << "] = " << arrayZ[i];
      }
      axisInPixelZ = new TAxis(3, arrayZ);

      hHitMap = new TH2F(Form("hHitMap_R%.0f", layerRadius), "HitMap;z (cm); r#phi (cm)", 1000, -1000, 1000, 1000, -1000, 1000);
      hHitMapInPixel = new TH2F(Form("hHitMapInPixel_R%.0f", layerRadius), "HitMapInPixel;z (cm); r#phi (cm)", 1000, -10, 10, 1000, -10, 10);
      hHitMapInPixelBefore = new TH2F(Form("hHitMapInPixelBefore_R%.0f", layerRadius), "HitMapInPixel;z (cm); r#phi (cm)", 1000, -10, 10, 1000, -10, 10);
    }

    bool isInTOFActiveArea(std::array<float, 3> hitPosition)
    {
      if (fractionInactive <= 0.0f) {
        return true;
      }
      if (fractionInactive >= 1.0f) {
        return false;
      }

      // Convert 3D position to cylindrical coordinates for area calculation
      const float phi = std::atan2(hitPosition[1], hitPosition[0]);
      const float rphi = phi * layerRadius;
      const float z = hitPosition[2];
      const float r = std::sqrt(hitPosition[0] * hitPosition[0] + hitPosition[1] * hitPosition[1]);

      // Check if hit is within layer geometric acceptance
      if (std::abs(layerRadius - r) > 10.f) {
        LOG(debug) << "Hit out of TOF layer acceptance: r=" << r << " cm with respect to the layer radius " << layerRadius;
        return false;
      }
      if (std::abs(z) > layerLength / 2.0f) {
        LOG(debug) << "Hit out of TOF layer acceptance: z=" << z << " cm with respect to the layer length " << layerLength;
        return false;
      }

      const int pixelIndexPhi = axisRPhi->FindBin(rphi);
      const int pixelIndexZ = axisZ->FindBin(z);

      // LOG(info) << "Hit pixel " << pixelIndexPhi << "/" << nPixelsRPhi << " and " << pixelIndexZ << "/" << nPixelsZ;

      if (pixelIndexPhi <= 0 || pixelIndexPhi > axisRPhi->GetNbins() || pixelIndexZ <= 0 || pixelIndexZ > axisZ->GetNbins()) {
        // LOG(info) << "Hit out of TOF layer pixel range: pixelIndexPhi=" << pixelIndexPhi << ", pixelIndexZ=" << pixelIndexZ;
        return false;
      }
      // Calculate local position within the pixel
      const float localRPhi = (rphi - axisRPhi->GetBinCenter(pixelIndexPhi));
      const float localZ = (z - axisZ->GetBinCenter(pixelIndexZ));

      // The difference between the hit and the pixel position cannot be greater than the size of the pixel
      if (std::abs(localRPhi - axisRPhi->GetBinCenter(pixelIndexPhi)) > axisRPhi->GetBinWidth(pixelIndexPhi)) {
        // LOG(warning) << "Local hit difference in phi is bigger than the pixel size";
      }
      if (std::abs(localZ - axisZ->GetBinCenter(pixelIndexZ)) > axisZ->GetBinWidth(pixelIndexZ)) {
        // LOG(warning) << "Local hit difference in z is bigger than the pixel size";
      }
      hHitMapInPixelBefore->Fill(localZ, localRPhi);
      switch (axisInPixelRPhi->FindBin(localRPhi)) {
        case 0:
        case 1:
        case 3:
        case 4:
          return false;
      }
      switch (axisInPixelZ->FindBin(localZ)) {
        case 0:
        case 1:
        case 3:
        case 4:
          return false;
      }
      hHitMapInPixel->Fill(localZ, localRPhi);
      hHitMap->Fill(z, rphi);
      return true;
    }

    /// Check if a track hits the active area (convenience function)
    /// \param track the track parameters for automatic hit position calculation
    /// \return true if the hit is in the active area
    bool isTrackInActiveArea(const o2::track::TrackParCov& track)
    {
      if (fractionInactive <= 0.f) {
        return true;
      }
      float x, y, z;
      if (!track.getXatLabR(layerRadius, x, magField)) {
        LOG(debug) << "Could not propagate track to TOF layer at radius " << layerRadius << " cm";
        return false;
      }
      bool b;
      ROOT::Math::PositionVector3D hit = track.getXYZGloAt(x, magField, b);
      if (!b) {
        LOG(debug) << "Could not get hit position at radius " << layerRadius << " cm";
        return false;
      }
      hit.GetCoordinates(x, y, z);
      return isInTOFActiveArea(std::array<float, 3>{x, y, z});
    }
  };

  bool isInInnerTOFActiveArea(const o2::track::TrackParCov& track)
  {
    static TOFLayerEfficiency innerTOFLayer(simConfig.innerTOFRadius, simConfig.innerTOFLength, simConfig.innerTOFPixelDimension, simConfig.innerTOFFractionOfInactiveArea, mMagneticField);
    return innerTOFLayer.isTrackInActiveArea(track);
  }

  bool isInOuterTOFActiveArea(const o2::track::TrackParCov& track)
  {
    static TOFLayerEfficiency outerTOFLayer(simConfig.outerTOFRadius, simConfig.outerTOFLength, simConfig.outerTOFPixelDimension, simConfig.outerTOFFractionOfInactiveArea, mMagneticField);
    return outerTOFLayer.isTrackInActiveArea(track);
  }

  /// function to calculate track length of this track up to a certain radius
  /// \param track the input track
  /// \param radius the radius of the layer you're calculating the length to
  /// \param magneticField the magnetic field to use when propagating
  static float computeTrackLength(o2::track::TrackParCov track, float radius, float magneticField)
  {
    // don't make use of the track parametrization
    float length = -100;

    o2::math_utils::CircleXYf_t trcCircle;
    float sna, csa;
    track.getCircleParams(magneticField, trcCircle, sna, csa);

    // distance between circle centers (one circle is at origin -> easy)
    const float centerDistance = std::hypot(trcCircle.xC, trcCircle.yC);

    // condition of circles touching - if not satisfied returned length will be -100
    if (centerDistance < trcCircle.rC + radius && centerDistance > std::fabs(trcCircle.rC - radius)) {
      length = 0.0f;

      // base radical direction
      const float ux = trcCircle.xC / centerDistance;
      const float uy = trcCircle.yC / centerDistance;
      // calculate perpendicular vector (normalized) for +/- displacement
      const float vx = -uy;
      const float vy = +ux;
      // calculate coordinate for radical line
      const float radical = (centerDistance * centerDistance - trcCircle.rC * trcCircle.rC + radius * radius) / (2.0f * centerDistance);
      // calculate absolute displacement from center-to-center axis
      const float displace = (0.5f / centerDistance) * std::sqrt(
                                                         (-centerDistance + trcCircle.rC - radius) *
                                                         (-centerDistance - trcCircle.rC + radius) *
                                                         (-centerDistance + trcCircle.rC + radius) *
                                                         (centerDistance + trcCircle.rC + radius));

      // possible intercept points of track and TOF layer in 2D plane
      const float point1[2] = {radical * ux + displace * vx, radical * uy + displace * vy};
      const float point2[2] = {radical * ux - displace * vx, radical * uy - displace * vy};

      // decide on correct intercept point
      std::array<float, 3> mom;
      track.getPxPyPzGlo(mom);
      const float scalarProduct1 = point1[0] * mom[0] + point1[1] * mom[1];
      const float scalarProduct2 = point2[0] * mom[0] + point2[1] * mom[1];

      // get start point
      std::array<float, 3> startPoint;
      track.getXYZGlo(startPoint);

      float cosAngle = -1000, modulus = -1000;

      if (scalarProduct1 > scalarProduct2) {
        modulus = std::hypot(point1[0] - trcCircle.xC, point1[1] - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
        cosAngle = (point1[0] - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (point1[1] - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
      } else {
        modulus = std::hypot(point2[0] - trcCircle.xC, point2[1] - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
        cosAngle = (point2[0] - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (point2[1] - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
      }
      cosAngle /= modulus;
      length = trcCircle.rC * std::acos(cosAngle);
      length *= std::sqrt(1.0f + track.getTgl() * track.getTgl());
    }
    return length;
  }

  /// returns velocity in centimeters per picoseconds
  /// \param momentum the momentum of the tarck
  /// \param mass the mass of the particle
  float computeParticleVelocity(float momentum, float mass)
  {
    const float a = momentum / mass;
    // uses light speed in cm/ps so output is in those units
    return o2::constants::physics::LightSpeedCm2PS * a / std::sqrt((1.f + a * a));
  }

  struct TracksWithTime {
    TracksWithTime(int pdgCode,
                   std::pair<float, float> innerTOFTime,
                   std::pair<float, float> outerTOFTime,
                   std::pair<float, float> trackLengthInnerTOF,
                   std::pair<float, float> trackLengthOuterTOF,
                   std::pair<float, float> momentum,
                   std::pair<float, float> pseudorapidity,
                   float noSmearingPt)
      : mPdgCode(pdgCode),
        mInnerTOFTime(innerTOFTime),
        mOuterTOFTime(outerTOFTime),
        mTrackLengthInnerTOF(trackLengthInnerTOF),
        mTrackLengthOuterTOF(trackLengthOuterTOF),
        mMomentum(momentum),
        mPseudorapidity{pseudorapidity},
        mNoSmearingPt{noSmearingPt} {}

    int mPdgCode;
    std::pair<float, float> mInnerTOFTime;        // Measured time and expected resolution at the inner TOF [ps]
    std::pair<float, float> mOuterTOFTime;        // Measured time and expected resolution at the outer TOF [ps]
    std::pair<float, float> mTrackLengthInnerTOF; // Measured and generated length at the innerTOF [cm]
    std::pair<float, float> mTrackLengthOuterTOF; // Measured and generated length at the outerTOF [cm]
    std::pair<float, float> mMomentum;            // Measured momentum and uncertainty on the momentum [GeV/c]
    std::pair<float, float> mPseudorapidity;      // Measured pseudorapidity and uncertainty on the pseudorapidity
    float mNoSmearingPt;                          // No smearing pt
  };

  std::vector<TracksWithTime> tracksWithTime;
  bool eventTime(std::vector<TracksWithTime>& tracks,
                 std::array<float, 2>& tzero)
  {

    float sum = 0.;
    float sumw = 0.;

    // Todo: check the different mass hypothesis iteratively
    for (const auto& track : tracks) {
      auto pdgInfo = pdg->GetParticle(track.mPdgCode);
      if (pdgInfo == nullptr) {
        continue;
      }
      const float mass = pdgInfo->Mass();
      const float mass2 = mass * mass;
      const float tof = track.mInnerTOFTime.first; // [ps]
      if (tof <= 0) {
        continue;
      }
      const float etof = track.mInnerTOFTime.second;         // [ps]
      const float length = track.mTrackLengthInnerTOF.first; // [cm]
      float p = track.mMomentum.first;                       // [GeV/c]
      p *= std::abs(pdgInfo->Charge()) / 3.;                 // Total momentum
      const float ep = track.mMomentum.second;               // [GeV/c]
      const float p2 = p * p;
      const float lengthOverC = length * o2::constants::physics::invLightSpeedCm2PS;
      const float texp = lengthOverC / p * std::sqrt(mass2 + p2);
      // LOG(info) << "TOF: " << tof << " " << etof << " vs " << texp;
      const float etexp = lengthOverC * mass2 / p2 / std::sqrt(mass2 + p2) * ep;
      const float sigma = std::sqrt(etexp * etexp + etof * etof);
      const float deltat = tof - texp;

      const float w = 1. / (sigma * sigma);

      sum += w * deltat;
      sumw += w;
    }

    static constexpr float kMaxEventTimeResolution = 200.f;
    if (sumw <= 0. || tracks.size() <= 1 || std::sqrt(1. / sumw) > kMaxEventTimeResolution) {
      tzero[0] = 0.;                      // [ps]
      tzero[1] = kMaxEventTimeResolution; // [ps]
      return false;
    }

    tzero[0] = sum / sumw;
    tzero[1] = std::sqrt(1. / sumw);
    return true;
  }

  /// returns track time resolution
  /// \param pt the transverse momentum of the tarck
  /// \param eta the pseudorapidity of the tarck
  /// \param trackPtResolution the absolute resolution on pt
  /// \param trackEtaResolution the absolute resolution on eta
  /// \param mass the mass of the particle
  /// \param detRadius the radius of the cylindrical layer
  /// \param magneticField the magnetic field (along Z)
  double calculateTrackTimeResolutionAdvanced(const float pt,
                                              const float eta,
                                              const float trackPtResolution,
                                              const float trackEtaResolution,
                                              const float mass,
                                              const float detRadius,
                                              const float magneticField)
  {
    // Compute tracking contribution to timing using the error propagation formula
    // Uses light speed in m/ps, magnetic field in T (*0.1 for conversion kGauss -> T)
    double a0 = mass * mass;
    double a1 = 0.299792458 * (0.1 * magneticField) * (0.01 * o2::constants::physics::LightSpeedCm2NS / 1e+3);
    double a2 = (detRadius * 0.01) * (detRadius * 0.01) * (0.299792458) * (0.299792458) * (0.1 * magneticField) * (0.1 * magneticField) / 2.0;
    double dtofOndPt = (std::pow(pt, 4) * std::pow(std::cosh(eta), 2) * std::acos(1.0 - a2 / std::pow(pt, 2)) - 2.0 * a2 * std::pow(pt, 2) * (a0 + std::pow(pt * std::cosh(eta), 2)) / std::sqrt(a2 * (2.0 * std::pow(pt, 2) - a2))) / (a1 * std::pow(pt, 3) * std::sqrt(a0 + std::pow(pt * std::cosh(eta), 2)));
    double dtofOndEta = std::pow(pt, 2) * std::sinh(eta) * std::cosh(eta) * std::acos(1.0 - a2 / std::pow(pt, 2)) / (a1 * std::sqrt(a0 + std::pow(pt * std::cosh(eta), 2)));
    double trackTimeResolution = std::hypot(std::fabs(dtofOndPt) * trackPtResolution, std::fabs(dtofOndEta) * trackEtaResolution);
    return trackTimeResolution;
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::OTFLUTConfigId>::iterator const& collision,
               soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels> const& tracks,
               aod::McParticles const&,
               aod::McCollisions const&)
  {
    o2::dataformats::VertexBase pvVtx({collision.posX(), collision.posY(), collision.posZ()},
                                      {collision.covXX(), collision.covXY(), collision.covYY(),
                                       collision.covXZ(), collision.covYZ(), collision.covZZ()});

    std::array<float, 6> mcPvCov = {0.};
    o2::dataformats::VertexBase mcPvVtx({0.0f, 0.0f, 0.0f}, mcPvCov);
    const float eventCollisionTimePS = (simConfig.considerEventTime.value ? collision.collisionTime() * 1e3 : 0.f); // convert ns to ps
    if (collision.has_mcCollision()) {
      auto mcCollision = collision.mcCollision();
      mcPvVtx.setX(mcCollision.posX());
      mcPvVtx.setY(mcCollision.posY());
      mcPvVtx.setZ(mcCollision.posZ());
    } // else remains untreated for now

    // First we compute the number of charged particles in the event if LUTs are loaded
    float dNdEta = 0.f;
    for (const auto& track : tracks) {
      if (!track.has_mcParticle())
        continue;
      auto mcParticle = track.mcParticle();
      if (std::abs(mcParticle.eta()) > simConfig.multiplicityEtaRange) {
        continue;
      }
      if (mcParticle.has_daughters()) {
        continue;
      }
      const auto& pdgInfo = pdg->GetParticle(mcParticle.pdgCode());
      if (!pdgInfo) {
        // LOG(warning) << "PDG code " << mcParticle.pdgCode() << " not found in the database";
        continue;
      }
      if (pdgInfo->Charge() == 0) {
        continue;
      }
      dNdEta += 1.f;
    }
    if (plotsConfig.doQAplots) {
      histos.fill(HIST("h1dNdeta"), dNdEta);
    }

    tracksWithTime.clear(); // clear the vector of tracks with time to prepare the cache for the next event
    tracksWithTime.reserve(tracks.size());
    // First loop to generate the arrival time of the particles to the TOF layers
    for (const auto& track : tracks) {
      // first step: find precise arrival time (if any)
      // --- convert track into perfect track
      if (!track.has_mcParticle()) { // should always be OK but check please
        upgradeTofMC(-999.f, -999.f, -999.f, -999.f);
        continue;
      } else {
        LOG(debug) << "Track without mcParticle found!";
      }
      const auto& mcParticle = track.mcParticle();
      o2::track::TrackParCov o2track = o2::upgrade::convertMCParticleToO2Track(mcParticle, pdg);

      float xPv = -100.f;
      static constexpr float kTrkXThreshold = -99.f; // Threshold to consider a good propagation of the track
      if (o2track.propagateToDCA(mcPvVtx, mMagneticField)) {
        xPv = o2track.getX();
      }
      float trackLengthInnerTOF = -1, trackLengthOuterTOF = -1;
      if (xPv > kTrkXThreshold) {
        trackLengthInnerTOF = computeTrackLength(o2track, simConfig.innerTOFRadius, mMagneticField);
        trackLengthOuterTOF = computeTrackLength(o2track, simConfig.outerTOFRadius, mMagneticField);
      }

      // Check if the track hit a sensitive area of the TOF
      const bool activeInnerTOF = isInInnerTOFActiveArea(o2track);
      if (!activeInnerTOF) {
        trackLengthInnerTOF = -999.f;
      } else {
        float x = 0.f;
        o2track.getXatLabR(simConfig.innerTOFRadius, x, mMagneticField);
        histos.fill(HIST("iTOF/h2HitMap"), o2track.getZAt(x, mMagneticField), simConfig.innerTOFRadius * o2track.getPhiAt(x, mMagneticField));
      }

      const bool activeOuterTOF = isInOuterTOFActiveArea(o2track);
      if (!activeOuterTOF) {
        trackLengthOuterTOF = -999.f;
      } else {
        float x = 0.f;
        o2track.getXatLabR(simConfig.outerTOFRadius, x, mMagneticField);
        histos.fill(HIST("oTOF/h2HitMap"), o2track.getZAt(x, mMagneticField), simConfig.outerTOFRadius * o2track.getPhiAt(x, mMagneticField));
      }

      // get mass to calculate velocity
      auto pdgInfo = pdg->GetParticle(mcParticle.pdgCode());
      if (pdgInfo == nullptr) {
        LOG(error) << "PDG code " << mcParticle.pdgCode() << " not found in the database";
        upgradeTofMC(-999.f, -999.f, -999.f, -999.f);
        continue;
      }
      const float v = computeParticleVelocity(o2track.getP(), pdgInfo->Mass());
      const float expectedTimeInnerTOF = trackLengthInnerTOF > 0 ? trackLengthInnerTOF / v + eventCollisionTimePS : -999.f; // arrival time to the Inner TOF in ps
      const float expectedTimeOuterTOF = trackLengthOuterTOF > 0 ? trackLengthOuterTOF / v + eventCollisionTimePS : -999.f; // arrival time to the Outer TOF in ps
      upgradeTofMC(expectedTimeInnerTOF, trackLengthInnerTOF, expectedTimeOuterTOF, trackLengthOuterTOF);

      // Smear with expected resolutions
      const float measuredTimeInnerTOF = pRandomNumberGenerator.Gaus(expectedTimeInnerTOF, simConfig.innerTOFTimeReso);
      const float measuredTimeOuterTOF = pRandomNumberGenerator.Gaus(expectedTimeOuterTOF, simConfig.outerTOFTimeReso);

      // Now we calculate the expected arrival time following certain mass hypotheses
      // and the (imperfect!) reconstructed track parametrizations
      float trackLengthRecoInnerTOF = -1, trackLengthRecoOuterTOF = -1;
      auto recoTrack = getTrackParCov(track);
      xPv = -100.f;
      if (recoTrack.propagateToDCA(pvVtx, mMagneticField)) {
        xPv = recoTrack.getX();
      }
      if (xPv > kTrkXThreshold) {
        trackLengthRecoInnerTOF = computeTrackLength(recoTrack, simConfig.innerTOFRadius, mMagneticField);
        trackLengthRecoOuterTOF = computeTrackLength(recoTrack, simConfig.outerTOFRadius, mMagneticField);
      }

      // cache the track info needed for the event time calculation
      // Reuse or emplace a new object in the vector
      tracksWithTime.emplace_back(TracksWithTime{mcParticle.pdgCode(),
                                                 {measuredTimeInnerTOF, simConfig.innerTOFTimeReso},
                                                 {measuredTimeOuterTOF, simConfig.outerTOFTimeReso},
                                                 {trackLengthRecoInnerTOF, trackLengthInnerTOF},
                                                 {trackLengthRecoOuterTOF, trackLengthOuterTOF},
                                                 {recoTrack.getP(), recoTrack.getSigma1Pt2()},
                                                 {recoTrack.getEta(), recoTrack.getSigmaTgl2()},
                                                 o2track.getPt()});
    }

    // Now we compute the event time for the tracks

    std::array<float, 2> tzero = {0.f, 0.f};
    bool etStatus = false;
    if (simConfig.considerEventTime.value) {
      etStatus = eventTime(tracksWithTime, tzero);
      if (!etStatus) {
        LOG(debug) << "Event time calculation failed with " << tracksWithTime.size() << " tracks with time and " << dNdEta << " charged particles";
      }
    }

    if (plotsConfig.doQAplots) {
      histos.fill(HIST("h2dEventTime"), tzero[0], eventCollisionTimePS);
      histos.fill(HIST("h1dEventTimegen"), eventCollisionTimePS);
      histos.fill(HIST("h1dEventTimerec"), tzero[0]);
      histos.fill(HIST("h2dEventTimeres"), dNdEta, tzero[1]);
      if (etStatus) {
        histos.fill(HIST("h1dEventTimedelta"), eventCollisionTimePS - tzero[0]);
      }
      static_cast<TEfficiency*>(listEfficiency->At(0))->Fill(etStatus, dNdEta);
    }

    // Then we do a second loop to compute the measured quantities with the measured event time
    int trackWithTimeIndex = 0;
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) { // should always be OK but check please
        continue;
      }
      const auto& mcParticle = track.mcParticle();

      const auto& trkWithTime = tracksWithTime[trackWithTimeIndex++];
      const float trackLengthRecoInnerTOF = trkWithTime.mTrackLengthInnerTOF.first;
      const float trackLengthRecoOuterTOF = trkWithTime.mTrackLengthOuterTOF.first;
      const float trackLengthInnerTOF = trkWithTime.mTrackLengthInnerTOF.second;
      const float trackLengthOuterTOF = trkWithTime.mTrackLengthOuterTOF.second;
      // Todo: remove the bias of the track used in the event time calculation for low multiplicity events
      const float measuredTimeInnerTOF = trkWithTime.mInnerTOFTime.first - tzero[0];
      const float measuredTimeOuterTOF = trkWithTime.mOuterTOFTime.first - tzero[0];
      const float momentum = trkWithTime.mMomentum.first;
      const float pseudorapidity = trkWithTime.mPseudorapidity.first;
      const float noSmearingPt = trkWithTime.mNoSmearingPt;

      // Straight to Nsigma
      static std::array<float, kParticles> expectedTimeInnerTOF, expectedTimeOuterTOF;
      static std::array<float, kParticles> deltaTimeInnerTOF, deltaTimeOuterTOF;
      static std::array<float, kParticles> nSigmaInnerTOF, nSigmaOuterTOF;
      static constexpr int kParticlePdgs[kParticles] = {kElectron,
                                                        kMuonMinus,
                                                        kPiPlus,
                                                        kKPlus,
                                                        kProton,
                                                        o2::constants::physics::kDeuteron,
                                                        o2::constants::physics::kTriton,
                                                        o2::constants::physics::kHelium3,
                                                        o2::constants::physics::kAlpha};
      static constexpr float kParticleMasses[kParticles] = {o2::constants::physics::MassElectron,
                                                            o2::constants::physics::MassMuon,
                                                            o2::constants::physics::MassPionCharged,
                                                            o2::constants::physics::MassKaonCharged,
                                                            o2::constants::physics::MassProton,
                                                            o2::constants::physics::MassDeuteron,
                                                            o2::constants::physics::MassTriton,
                                                            o2::constants::physics::MassHelium3,
                                                            o2::constants::physics::MassAlpha};
      static constexpr float kParticleCharges[kParticles] = {1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 1.f, 2.f, 2.f};
      float momentumHypotheses[kParticles]; // Store momentum hypothesis for each particle

      auto truePdgInfo = pdg->GetParticle(mcParticle.pdgCode());
      float rigidity = momentum; // fallback to momentum if charge unknown

      // Use MC truth charge for rigidity calculation
      if (truePdgInfo) {
        const float trueCharge = std::abs(truePdgInfo->Charge()) / 3.0f;
        if (trueCharge > 0) {
          rigidity = momentum / trueCharge;
        }
      }

      if (plotsConfig.doQAplots) {
        // unit conversion: length in cm, time in ps
        const float innerBeta = (trackLengthInnerTOF / measuredTimeInnerTOF) / o2::constants::physics::LightSpeedCm2PS;
        const float outerBeta = (trackLengthOuterTOF / measuredTimeOuterTOF) / o2::constants::physics::LightSpeedCm2PS;
        if (trackLengthRecoInnerTOF > 0) {
          histos.fill(HIST("iTOF/h2dVelocityVsMomentumInner"), momentum, innerBeta);
          histos.fill(HIST("iTOF/h2dVelocityVsRigidityInner"), rigidity, innerBeta);
          histos.fill(HIST("iTOF/h2dTrackLengthInnerVsPt"), noSmearingPt, trackLengthInnerTOF);
          histos.fill(HIST("iTOF/h2dTrackLengthInnerRecoVsPt"), noSmearingPt, trackLengthRecoInnerTOF);
        }
        if (trackLengthRecoOuterTOF > 0) {
          histos.fill(HIST("oTOF/h2dVelocityVsMomentumOuter"), momentum, outerBeta);
          histos.fill(HIST("oTOF/h2dVelocityVsRigidityOuter"), rigidity, outerBeta);
          histos.fill(HIST("oTOF/h2dTrackLengthOuterVsPt"), noSmearingPt, trackLengthOuterTOF);
          histos.fill(HIST("oTOF/h2dTrackLengthOuterRecoVsPt"), noSmearingPt, trackLengthRecoOuterTOF);
        }
      }

      // For every mass hypothesis compute the expected time, the delta with respect to it and the nsigma
      for (int ii = 0; ii < kParticles; ii++) {
        expectedTimeInnerTOF[ii] = -100;
        expectedTimeOuterTOF[ii] = -100;
        deltaTimeInnerTOF[ii] = -100;
        deltaTimeOuterTOF[ii] = -100;
        nSigmaInnerTOF[ii] = -100;
        nSigmaOuterTOF[ii] = -100;

        momentumHypotheses[ii] = rigidity * kParticleCharges[ii]; // Total momentum for this hypothesis
        const float v = computeParticleVelocity(momentumHypotheses[ii], kParticleMasses[ii]);

        expectedTimeInnerTOF[ii] = trackLengthInnerTOF / v;
        expectedTimeOuterTOF[ii] = trackLengthOuterTOF / v;

        deltaTimeInnerTOF[ii] = measuredTimeInnerTOF - expectedTimeInnerTOF[ii];
        deltaTimeOuterTOF[ii] = measuredTimeOuterTOF - expectedTimeOuterTOF[ii];

        // Evaluate total sigma (layer + tracking resolution)
        float innerTotalTimeReso = simConfig.innerTOFTimeReso;
        float outerTotalTimeReso = simConfig.outerTOFTimeReso;
        if (simConfig.flagIncludeTrackTimeRes) {
          const float transverseMomentum = momentumHypotheses[ii] / std::cosh(pseudorapidity);
          double ptResolution = transverseMomentum * transverseMomentum * std::sqrt(trkWithTime.mMomentum.second);
          double etaResolution = std::fabs(std::sin(2.0 * std::atan(std::exp(-pseudorapidity)))) * std::sqrt(trkWithTime.mPseudorapidity.second);
          if (simConfig.flagTOFLoadDelphesLUTs) {
            if (mSmearer[collision.lutConfigId()]->hasTable(kParticlePdgs[ii])) { // Only if the LUT for this particle was loaded
              ptResolution = mSmearer[collision.lutConfigId()]->getAbsPtRes(kParticlePdgs[ii], dNdEta, pseudorapidity, transverseMomentum);
              etaResolution = mSmearer[collision.lutConfigId()]->getAbsEtaRes(kParticlePdgs[ii], dNdEta, pseudorapidity, transverseMomentum);
            }
          }
          const float innerTrackTimeReso = calculateTrackTimeResolutionAdvanced(transverseMomentum, pseudorapidity, ptResolution, etaResolution, kParticleMasses[ii], simConfig.innerTOFRadius, mMagneticField);
          const float outerTrackTimeReso = calculateTrackTimeResolutionAdvanced(transverseMomentum, pseudorapidity, ptResolution, etaResolution, kParticleMasses[ii], simConfig.outerTOFRadius, mMagneticField);
          innerTotalTimeReso = std::hypot(simConfig.innerTOFTimeReso, innerTrackTimeReso);
          outerTotalTimeReso = std::hypot(simConfig.outerTOFTimeReso, outerTrackTimeReso);

          if (plotsConfig.doQAplots) {
            if (std::fabs(mcParticle.pdgCode()) == kParticlePdgs[ii]) {
              if (trackLengthRecoInnerTOF > 0) {
                h2dInnerTimeResTrack[ii]->Fill(momentumHypotheses[ii], innerTrackTimeReso);
                h2dInnerTimeResTotal[ii]->Fill(momentumHypotheses[ii], innerTotalTimeReso);
              }
              if (trackLengthRecoOuterTOF > 0) {
                h2dOuterTimeResTrack[ii]->Fill(momentumHypotheses[ii], outerTrackTimeReso);
                h2dOuterTimeResTotal[ii]->Fill(momentumHypotheses[ii], outerTotalTimeReso);
                static constexpr int kIdPion = 2;
                if (ii == kIdPion) {
                  histos.fill(HIST("h2dRelativePtResolution"), transverseMomentum, 100.0 * ptResolution / transverseMomentum);
                  histos.fill(HIST("h2dRelativeEtaResolution"), pseudorapidity, 100.0 * etaResolution / (std::fabs(pseudorapidity) + 1e-6));
                }
              }
            }
          }
        }

        // Fixme: assumes dominant resolution effect is the TOF resolution
        // and not the tracking itself. It's *probably* a fair assumption
        // but it should be tested further! --> FIXED IN THIS VERSION
        if (trackLengthInnerTOF > 0 && trackLengthRecoInnerTOF > 0)
          nSigmaInnerTOF[ii] = deltaTimeInnerTOF[ii] / std::sqrt(innerTotalTimeReso * innerTotalTimeReso + tzero[1] * tzero[1]);
        if (trackLengthOuterTOF > 0 && trackLengthRecoOuterTOF > 0)
          nSigmaOuterTOF[ii] = deltaTimeOuterTOF[ii] / std::sqrt(outerTotalTimeReso * outerTotalTimeReso + tzero[1] * tzero[1]);
      }

      if (plotsConfig.doQAplots) {
        for (int ii = 0; ii < kParticles; ii++) {
          if (std::fabs(mcParticle.pdgCode()) != pdg->GetParticle(kParticlePdgs[ii])->PdgCode()) {
            continue;
          }
          if (trackLengthRecoInnerTOF > 0) {
            for (int iii = 0; iii < kParticles; iii++) {
              h2dInnerNsigmaTrue[ii][iii]->Fill(momentumHypotheses[ii], nSigmaInnerTOF[iii]);
              h2dInnerDeltaTrue[ii][iii]->Fill(momentumHypotheses[ii], deltaTimeInnerTOF[iii]);
            }
          }
          if (trackLengthRecoOuterTOF > 0) {
            for (int iii = 0; iii < kParticles; iii++) {
              h2dOuterNsigmaTrue[ii][iii]->Fill(momentumHypotheses[ii], nSigmaOuterTOF[iii]);
              h2dOuterDeltaTrue[ii][iii]->Fill(momentumHypotheses[ii], deltaTimeOuterTOF[iii]);
            }
          }
        }
      }

      const float deltaTrackLengthInnerTOF = std::abs(trackLengthInnerTOF - trackLengthRecoInnerTOF);
      if (trackLengthInnerTOF > 0 && trackLengthRecoInnerTOF > 0) {
        histos.fill(HIST("iTOF/h2dDeltaTrackLengthInnerVsPt"), noSmearingPt, deltaTrackLengthInnerTOF);
      }
      const float deltaTrackLengthOuterTOF = std::abs(trackLengthOuterTOF - trackLengthRecoOuterTOF);
      if (trackLengthOuterTOF > 0 && trackLengthRecoOuterTOF > 0) {
        histos.fill(HIST("oTOF/h2dDeltaTrackLengthOuterVsPt"), noSmearingPt, deltaTrackLengthOuterTOF);
      }

      // Sigmas have been fully calculated. Please populate the NSigma helper table (once per track)
      upgradeTof(tzero[0], tzero[1],
                 nSigmaInnerTOF[0], nSigmaInnerTOF[1], nSigmaInnerTOF[2], nSigmaInnerTOF[3], nSigmaInnerTOF[4], nSigmaInnerTOF[5], nSigmaInnerTOF[6], nSigmaInnerTOF[7], nSigmaInnerTOF[8],
                 measuredTimeInnerTOF, trackLengthRecoInnerTOF,
                 nSigmaOuterTOF[0], nSigmaOuterTOF[1], nSigmaOuterTOF[2], nSigmaOuterTOF[3], nSigmaOuterTOF[4], nSigmaOuterTOF[5], nSigmaOuterTOF[6], nSigmaOuterTOF[7], nSigmaOuterTOF[8],
                 measuredTimeOuterTOF, trackLengthRecoOuterTOF);
      upgradeTofExpectedTime(expectedTimeInnerTOF[0], expectedTimeInnerTOF[1], expectedTimeInnerTOF[2], expectedTimeInnerTOF[3], expectedTimeInnerTOF[4], expectedTimeInnerTOF[5], expectedTimeInnerTOF[6], expectedTimeInnerTOF[7], expectedTimeInnerTOF[8],
                             expectedTimeOuterTOF[0], expectedTimeOuterTOF[1], expectedTimeOuterTOF[2], expectedTimeOuterTOF[3], expectedTimeOuterTOF[4], expectedTimeOuterTOF[5], expectedTimeOuterTOF[6], expectedTimeOuterTOF[7], expectedTimeOuterTOF[8]);
    }

    if (trackWithTimeIndex != tracks.size()) {
      LOG(fatal) << "Track with time index " << trackWithTimeIndex << " does not match the number of tracks " << tracks.size();
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<OnTheFlyTofPid>(cfgc)}; }
